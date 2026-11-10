"""Name-to-CaSilico workflow. Standard-library Python; external tools stay separate."""
from __future__ import annotations

from runtime import storage_path, settings
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import html
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.parse import quote, urlencode, urljoin, unquote
from urllib.request import Request, urlopen
import zipfile

ROOT = Path(__file__).resolve().parent
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ARCHIVE = "https://data.cyverse.org/dav-anon/iplant/projects/Bioinformatics_Resource_Centers/VEuPathDB/AmoebaDB/"
COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


def log(message):
    print(message, flush=True)


def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]


def slug(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.") or "sample"


def write_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"不是 FASTA 文件：{path}")
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def write_fasta(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def fetch(url, destination=None):
    """Atomic downloads. Cached files are reused only with a matching SHA256 sidecar."""
    if destination:
        destination = Path(destination)
        stamp = Path(str(destination) + ".source.json")
        if destination.exists() and stamp.exists():
            meta = json.loads(stamp.read_text(encoding="utf-8"))
            if meta.get("url") == url and meta.get("sha256") == checksum(destination):
                return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        try:
            req = Request(url, headers={"User-Agent": "AmoebaCrisp/0.1 (research data client)"})
            with urlopen(req, timeout=120) as response:
                if destination is None:
                    return response.read()
                temp = destination.with_name(destination.name + ".part")
                with temp.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
                expected = response.headers.get("Content-Length")
                if expected and temp.stat().st_size != int(expected):
                    raise OSError("下载长度与服务器记录不一致")
                temp.replace(destination)
            write_json(stamp, {"url": url, "sha256": checksum(destination),
                               "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat()})
            return destination
        except Exception:
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def checksum(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def directory_links(url):
    parser = Links()
    parser.feed(fetch(url).decode("utf-8"))
    return parser.links


def ncbi_records(name, limit):
    reports, token = [], ""
    while True:
        query = {"page_size": 1000}
        if token:
            query["page_token"] = token
        result = json.loads(fetch(f"{NCBI}/genome/taxon/{quote(name, safe='')}/dataset_report?{urlencode(query)}"))
        reports.extend(result.get("reports", []))
        token = result.get("next_page_token")
        if not token:
            break
    current = [r for r in reports if r.get("assembly_info", {}).get("assembly_status", "current") == "current"]
    # RefSeq and GenBank paired assemblies represent the same biological assembly.
    by_pair = {}
    for record in current:
        key = tuple(sorted([record["accession"], record.get("paired_accession", record["accession"])]))
        if key not in by_pair or record["accession"].startswith("GCF_"):
            by_pair[key] = record
    records = sorted(by_pair.values(), key=lambda r: (not bool(r.get("annotation_info")), r["accession"]))
    if not records:
        raise RuntimeError(f"NCBI 未找到当前基因组：{name}")
    if limit and len(records) > limit:
        log(f"NCBI 找到 {len(records)} 个独立组装；本次按设置仅使用 {limit} 个。")
        records = records[:limit]
    return records


def ncbi_download(name, cache, limit):
    records = ncbi_records(name, limit)
    samples = []
    for r in records:
        accession = r["accession"]
        folder = cache / "ncbi" / accession
        types = ["GENOME_FASTA", "GENOME_GFF", "PROT_FASTA", "CDS_FASTA"]
        params = urlencode([("include_annotation_type", value) for value in types])
        url = f"{NCBI}/genome/accession/{accession}/download?{params}"
        log(f"下载 NCBI：{accession} / {r['organism']['organism_name']}")
        bundle = fetch(url, folder / "dataset.zip")
        with zipfile.ZipFile(bundle) as archive:
            # Copy only expected files; never extract arbitrary archive paths.
            names = archive.namelist()
            genome = [n for n in names if n.endswith("_genomic.fna") and "cds_from" not in n and "rna_from" not in n]
            annotation = [n for n in names if n.endswith("/genomic.gff")]
            if len(genome) != 1:
                raise RuntimeError(f"NCBI 包中的基因组文件不明确：{accession}")
            for member, target in [(genome[0], "genome.fna")] + ([(annotation[0], "annotation.gff3")] if annotation else []):
                with archive.open(member) as source, (folder / target).open("wb") as output:
                    shutil.copyfileobj(source, output)
        sample = {"id": accession, "name": r["organism"]["organism_name"],
                  "genome": str(folder / "genome.fna"),
                  "gff": str(folder / "annotation.gff3") if annotation else None,
                  "source": "ncbi", "accession": accession, "metadata": r}
        samples.append(sample)
    return samples


def amoebadb_download(name, cache, limit, base=ARCHIVE):
    """Public VEuPathDB archive; release is preserved, never labelled current."""
    words = name.strip().split()
    genus = words[0]
    index = urljoin(base.rstrip("/") + "/", quote(genus) + "/")
    links = directory_links(index)
    prefix = genus[0] + words[1] if len(words) > 1 else genus[0]
    matches = {}
    for link in links:
        filename = unquote(link.split("/")[-1])
        match = re.fullmatch(r"AmoebaDB-(\d+)_(.+)_Genome\.fasta(?:\.gz)?", filename)
        if not match:
            continue
        release, code = match.groups()
        if not code.lower().startswith(prefix.lower()):
            continue
        if len(words) > 2 and slug(" ".join(words[2:])).replace("_", "").lower() not in code.replace("-", "").lower():
            continue
        if code not in matches or int(release) > matches[code][0]:
            matches[code] = (int(release), link)
    if not matches:
        raise RuntimeError(f"AmoebaDB 公开存档未找到：{name}。可选择 NCBI；不会静默换库。")
    samples = []
    for code, (release, genome_link) in sorted(matches.items())[:limit or None]:
        folder = cache / "amoebadb" / f"AmoebaDB-{release}_{slug(code)}"
        stem = f"AmoebaDB-{release}_{code}"
        annotation = next((v for v in links if unquote(v.split("/")[-1]) in [stem + ".gff", stem + ".gff3", stem + ".gff.gz"]), None)
        log(f"下载 AmoebaDB 存档：{code}，release {release}（不保证是最新版本）")
        genome = fetch(urljoin(index, genome_link), folder / ("genome.fna.gz" if genome_link.endswith(".gz") else "genome.fna"))
        gff = fetch(urljoin(index, annotation), folder / ("annotation.gff3.gz" if annotation.endswith(".gz") else "annotation.gff3")) if annotation else None
        samples.append({"id": slug(stem), "name": name, "genome": str(genome),
                        "gff": str(gff) if gff else None, "source": "amoebadb-archive",
                        "release": release, "url": urljoin(index, genome_link)})
    return samples


def gff_attributes(text):
    result = {}
    for part in text.strip().split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            result[key] = unquote(value)
    return result


CODONS = dict(zip((a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def translate(seq):
    return "".join(CODONS.get(seq[i:i+3], "X") for i in range(0, len(seq) - 2, 3)).rstrip("*")


def extract_models(sample, outdir):
    """Translate spliced CDS for orthology, keep UNSPLICED genomic loci for DNA design."""
    genome = dict(read_fasta(sample["genome"]))
    if not genome:
        raise ValueError("空基因组")
    parents, attributes, cds = {}, {}, {}
    opener = gzip.open if str(sample["gff"]).endswith(".gz") else open
    with opener(sample["gff"], "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                raise ValueError("注释必须是 GFF3（九列）；不能直接使用 GTF")
            contig, _, kind, start, end, _, strand, phase, raw = fields
            attrs = gff_attributes(raw)
            if "ID" in attrs:
                parents[attrs["ID"]] = attrs.get("Parent", "").split(",")[0]
                attributes[attrs["ID"]] = attrs
            if kind != "CDS":
                continue
            for parent in attrs.get("Parent", attrs.get("ID", "")).split(","):
                if parent:
                    cds.setdefault(parent, []).append((contig, int(start), int(end), strand, phase, attrs))
    by_gene = {}
    for transcript, parts in cds.items():
        contigs = {p[0] for p in parts}
        strands = {p[3] for p in parts}
        if len(contigs) != 1 or len(strands) != 1 or next(iter(strands)) not in ("+", "-"):
            continue
        contig, strand = parts[0][0], parts[0][3]
        if contig not in genome:
            raise ValueError(f"GFF contig {contig} 不在基因组中")
        parts.sort(key=lambda p: p[1], reverse=strand == "-")
        gene = transcript
        seen = set()
        while parents.get(gene) and gene not in seen:
            seen.add(gene)
            gene = parents[gene]
        if any(p[5].get("pseudo", "false").lower() == "true" for p in parts):
            continue
        if any(p[1] < 1 or p[2] > len(genome[contig]) or p[1] > p[2] for p in parts):
            raise ValueError(f"GFF 坐标超出基因组：{transcript}")
        fragments = [genome[contig][p[1]-1:p[2]] for p in parts]
        if strand == "-":
            fragments = [rc(s) for s in fragments]
        coding = "".join(fragments)
        # Internal phases complete codons across exon boundaries; never trim each exon.
        phase = int(parts[0][4]) if parts[0][4].isdigit() else 0
        protein = translate(coding[phase:])
        if len(protein) < 30 or "*" in protein:
            continue
        start, end = min(p[1] for p in parts), max(p[2] for p in parts)
        locus = genome[contig][start-1:end]
        if strand == "-":
            locus = rc(locus)
        record = {"gene": gene, "transcript": transcript, "protein": protein,
                  "dna": locus, "contig": contig, "start": start, "end": end,
                  "strand": strand, "cds_length": len(coding), "cds_parts": len(parts)}
        if gene not in by_gene or len(protein) > len(by_gene[gene]["protein"]):
            by_gene[gene] = record
    models = {}
    for i, (_, record) in enumerate(sorted(by_gene.items()), 1):
        key = f"{sample['id']}__g{i:06d}"
        models[key] = record
    if not models:
        raise RuntimeError(f"{sample['id']} 未提取到可用编码基因，请核对 GFF3")
    outdir.mkdir(parents=True, exist_ok=True)
    write_fasta(outdir / (sample["id"] + ".faa"), ((key, r["protein"]) for key, r in models.items()))
    log(f"{sample['id']}：保留 {len(models)} 个基因的代表蛋白及连续基因组区域")
    return models


def linux_path(path):
    text = str(Path(path).resolve())
    match = re.match(r"^([A-Za-z]):[\\/](.*)$", text)
    return "/mnt/" + match[1].lower() + "/" + match[2].replace("\\", "/") if match else text


def external(command, logfile, cwd=None):
    log("运行：" + " ".join(str(x) for x in command))
    with open(logfile, "w", encoding="utf-8") as handle:
        result = subprocess.run([str(x) for x in command], cwd=cwd, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"外部工具失败（{result.returncode}），详见 {logfile}")


def linux_command(executable, args, config):
    envbin = config.get("linux_env_bin", executable.rsplit("/", 1)[0])
    storage = linux_path(config["storage"])
    command = ["env", "PATH=" + envbin + ":/usr/local/bin:/usr/bin:/bin",
               "APPTAINER_CACHEDIR=" + storage + "/software/container-cache",
               "APPTAINER_TMPDIR=" + storage + "/temp", "TMPDIR=" + storage + "/temp",
               executable] + list(args)
    return ["wsl", "-d", config.get("wsl_distribution", "Ubuntu"), "--"] + command if os.name == "nt" else command


def orthogroups(samples, models, run, config, threads):
    if len(samples) == 1:
        log("只有一个独立组装：跳过 OrthoFinder，结果仅代表该参考基因组。")
        return [(f"REF_{i:06d}", [key]) for i, key in enumerate(models, 1)]
    executable = config["orthofinder"]
    output = run / "orthofinder"
    args = ["-f", linux_path(run / "proteins"), "-t", str(threads), "-a", str(threads),
            "-og", "-o", linux_path(output)]
    if os.name == "nt":
        command = linux_command(executable, args, config)
    else:
        command = [executable, "-f", str(run / "proteins"), "-t", str(threads), "-a", str(threads), "-og", "-o", str(output)]
    external(command, run / "orthofinder.log")
    tables = list(output.rglob("Orthogroups.tsv"))
    if len(tables) != 1:
        raise RuntimeError("未找到唯一的 OrthoFinder Orthogroups.tsv")
    result = []
    with tables[0].open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            members = []
            for sample in samples:
                cell = row.get(sample["id"], "").strip()
                genes = [g.strip() for g in cell.split(",") if g.strip()]
                if len(genes) != 1 or genes[0] not in models:
                    break
                members.extend(genes)
            else:
                result.append((row["Orthogroup"], members))
    log(f"OrthoFinder：找到 {len(result)} 个覆盖所有样本的单拷贝同源群")
    return result


class WidgetParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.active, self.buff, self.widgets = False, [], []

    def handle_starttag(self, tag, attrs):
        if tag == "script" and dict(attrs).get("type") == "application/json":
            self.active, self.buff = True, []

    def handle_data(self, data):
        if self.active:
            self.buff.append(data)

    def handle_endtag(self, tag):
        if tag == "script" and self.active:
            self.widgets.append(json.loads("".join(self.buff)))
            self.active = False


def scalar(value):
    if isinstance(value, list):
        return scalar(value[0]) if len(value) == 1 else json.dumps(value, ensure_ascii=False)
    return value


def parse_casilico(path):
    parser = WidgetParser()
    parser.feed(Path(path).read_text(encoding="utf-8"))
    for widget in parser.widgets:
        table = widget.get("x", {})
        headers = [html.unescape(re.sub("<[^>]+>", "", h)).strip()
                   for h in re.findall(r"<th\b[^>]*>(.*?)</th>", table.get("container", ""), re.S)]
        if "Guide_RNA" not in headers:
            continue
        columns = table.get("data", [])
        # OffTarget=FALSE yields a single header row. Reject ambiguous layouts.
        if len(headers) != len(columns):
            raise ValueError(f"不支持的 CaSilico 列结构：{path}")
        if not columns:
            return []
        if len({len(col) for col in columns}) != 1:
            raise ValueError("CaSilico 结果列长度不一致")
        return [{headers[j] or "row": scalar(columns[j][i]) for j in range(len(headers))}
                for i in range(len(columns[0]))]
    # Empty CaSilico results may lack a table. Preserve as zero candidates.
    return []


def design(group, members, models, run, config, threads):
    folder = run / "casilico" / slug(group)
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / "targets.fasta"
    write_fasta(fasta, ((key, models[key]["dna"]) for key in members))
    script = folder / "design.R"
    rq = lambda x: json.dumps(str(x).replace("\\", "/"))
    script.write_text("\n".join([
        f".libPaths(c({rq(config['r_library'])}, .libPaths()))",
        "library(CaSilico)",
        f"setwd({rq(folder)})",
        "options(browser=function(...) invisible(NULL))",
        f"CaSilico(ResultsFolder='design', TargetFasta={rq(fasta)}, CrisprTypes='casV_A',",
        f" ConservationMethod=1, ConservationThreshold=1, OffTarget=FALSE, OffAsk=FALSE, Threads={threads})",
        "writeLines(as.character(packageVersion('CaSilico')), 'casilico_version.txt')",
        "writeLines('OK', 'completed.txt')",
    ]) + "\n", encoding="utf-8")
    external([config["rscript"], "--vanilla", script], folder / "design.log", cwd=folder)
    if not (folder / "completed.txt").exists():
        raise RuntimeError(f"CaSilico 未完成：{folder}")
    tables = list(folder.rglob("Spacers_Information.Html"))
    if len(tables) > 1:
        raise RuntimeError("CaSilico 输出表不唯一")
    return parse_casilico(tables[0]) if tables else []


def exact_sites(seq, guide, reverse=None):
    """TTTV + spacer on either strand. Coordinates are 1-based on genomic plus strand."""
    hits = []
    for strand, dna in [("+", seq), ("-", reverse if reverse is not None else rc(seq))]:
        pos = dna.find(guide)
        while pos != -1:
            pam = dna[max(0, pos-4):pos]
            if len(pam) == 4 and pam[:3] == "TTT" and pam[3] in "ACG":
                start, end = ((pos+1, pos+len(guide)) if strand == "+" else
                              (len(seq)-pos-len(guide)+1, len(seq)-pos))
                hits.append({"start": start, "end": end, "strand": strand, "pam": pam})
            pos = dna.find(guide, pos+1)
    return hits


def mismatch_sites(seq, guide, max_mismatches, reverse=None):
    """Exhaustive gapless Hamming search using k+1 disjoint exact seeds.

    At most k substitutions leave at least one seed intact. PAM is checked
    independently. Ambiguous reference bases never count as confirmed matches.
    """
    if max_mismatches == 0:
        return [{**hit, "mismatches": 0} for hit in exact_sites(seq, guide, reverse)]
    hits = []
    length = len(guide)
    boundaries = [i * length // (max_mismatches + 1) for i in range(max_mismatches + 2)]
    for strand, dna in [("+", seq), ("-", reverse if reverse is not None else rc(seq))]:
        seen = set()
        for left, right in zip(boundaries, boundaries[1:]):
            seed = guide[left:right]
            pos = dna.find(seed)
            while pos != -1:
                start = pos - left
                if 4 <= start <= len(dna) - length and start not in seen:
                    seen.add(start)
                    pam = dna[start-4:start]
                    if pam[:3] == "TTT" and pam[3] in "ACG":
                        target = dna[start:start+length]
                        if set(target) <= set("ACGT"):
                            mismatches = sum(a != b for a, b in zip(guide, target))
                            if mismatches <= max_mismatches:
                                a, b = ((start+1, start+length) if strand == "+" else
                                        (len(seq)-start-length+1, len(seq)-start))
                                hits.append({"start": a, "end": b, "strand": strand,
                                             "pam": pam, "mismatches": mismatches})
                pos = dna.find(seed, pos+1)
    return hits


def write_tsv(path, rows, fields):
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(candidates, samples, run, config, background, max_mismatches=3):
    guides = sorted({str(row["Guide_RNA"]).upper().replace("U", "T") for row in candidates})
    if any(not re.fullmatch("[ACGT]{20}", guide) for guide in guides):
        raise ValueError("CaSilico Cas12a 输出含非预期的 spacer 序列或长度")
    coverage, background_counts, coordinates = {g: set() for g in guides}, {g: 0 for g in guides}, []
    near_counts, background_hits = {g: 0 for g in guides}, []
    # Process each contig once in memory. Scan both strands for genuine PAM-bearing DNA sites.
    for sample in samples:
        log(f"核查完整基因组上的 spacer + TTTV：{sample['id']}")
        for contig, sequence in read_fasta(sample["genome"]):
            reverse = rc(sequence)
            for guide in guides:
                hits = exact_sites(sequence, guide, reverse)
                if hits:
                    coverage[guide].add(sample["id"])
                    coordinates.extend({"sample": sample["id"], "contig": contig, "spacer_dna": guide, **hit} for hit in hits)
    for filename in background:
        log(f"背景筛查：{filename}（0–{max_mismatches} 个碱基替换，TTTV，无插入缺失）")
        for contig, sequence in read_fasta(filename):
            reverse = rc(sequence)
            for guide in guides:
                hits = mismatch_sites(sequence, guide, max_mismatches, reverse)
                background_counts[guide] += sum(hit["mismatches"] == 0 for hit in hits)
                near_counts[guide] += len(hits)
                background_hits.extend({"background": str(filename), "contig": contig,
                                        "spacer_dna": guide, **hit} for hit in hits)
    by_guide = {}
    for candidate in candidates:
        guide = str(candidate["Guide_RNA"]).upper().replace("U", "T")
        by_guide.setdefault(guide, []).append(candidate)
    rows = []
    for number, guide in enumerate(guides, 1):
        raw = by_guide[guide]
        gc = 100 * sum(c in "GC" for c in guide) / len(guide)
        all_covered = len(coverage[guide]) == len(samples)
        structure_ok = any("Natural structure of DR" in str(r.get("Self_Complementarity_Free_Energy", ""))
                           and "Unnatural" not in str(r.get("Self_Complementarity_Free_Energy", "")) for r in raw)
        passes = all_covered and config["gc_min"] <= gc <= config["gc_max"] and structure_ok
        rows.append({"id": f"crRNA_{number:05d}", "spacer_dna_5to3": guide,
                     "spacer_rna_5to3": guide.replace("T", "U"),
                     "full_crrna_rna_5to3": config["direct_repeat_rna"] + guide.replace("T", "U"),
                     "direct_repeat_source": "CaSilico casV_A default; verify compatibility with chosen Cas12a",
                     "groups": ";".join(sorted({r["group"] for r in raw})),
                     "covered_genomes": len(coverage[guide]), "total_genomes": len(samples),
                     "coverage_percent": round(100 * len(coverage[guide])/len(samples), 2),
                     "gc_percent": round(gc, 2), "passes_target_filters": passes,
                     "background_exact_pam_hits": background_counts[guide] if background else "NOT_EVALUATED",
                     "background_near_pam_hits": near_counts[guide] if background else "NOT_EVALUATED",
                     "background_mismatch_screen": f"Hamming<={max_mismatches}; TTTV; no indels" if background else "NOT_EVALUATED",
                     "scope": "single_reference" if len(samples) == 1 else "sampled_assemblies",
                     "casilico_ids": ";".join(str(r.get("Guide_RNA_ID", "")) for r in raw)})
    rows.sort(key=lambda r: (-r["coverage_percent"], not r["passes_target_filters"], abs(r["gc_percent"] - 50)))
    fields = list(rows[0]) if rows else ["id", "spacer_rna_5to3", "full_crrna_rna_5to3", "coverage_percent", "background_exact_pam_hits"]
    write_tsv(run / "crRNA_designs.tsv", rows, fields)
    write_tsv(run / "target_sites.tsv", coordinates, ["sample", "contig", "spacer_dna", "start", "end", "strand", "pam"])
    write_tsv(run / "background_sites.tsv", background_hits,
              ["background", "contig", "spacer_dna", "start", "end", "strand", "pam", "mismatches"])
    selected = [r for r in rows if r["passes_target_filters"] and (not background or r["background_near_pam_hits"] == 0)]
    write_fasta(run / "crRNA_spacers.fasta", ((r["id"], r["spacer_rna_5to3"]) for r in selected))
    write_fasta(run / "crRNA_full_CaSilico_default.fasta", ((r["id"], r["full_crrna_rna_5to3"]) for r in selected))
    headings = ["编号", "spacer RNA（5′→3′）", "完整 RNA（CaSilico 默认 DR）", "覆盖率", "GC%", "背景近似匹配"]
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in [r["id"], r["spacer_rna_5to3"],
                   r["full_crrna_rna_5to3"], str(r["coverage_percent"])+"%", r["gc_percent"],
                   r["background_near_pam_hits"] if background else "未评估"]) + "</tr>" for r in selected)
    group_file = run / "groups.json"
    group_info = json.loads(group_file.read_text(encoding="utf-8")) if group_file.exists() else None
    scope = f"本次处理 {group_info['processed']} / {group_info['eligible']} 个符合区域长度条件的基因/同源群。" if group_info else ""
    report = f"""<!doctype html><meta charset='utf-8'><title>AmoebaCrisp 结果</title>
    <style>body{{font:16px system-ui;margin:40px;line-height:1.7}}td,th{{padding:8px;border:1px solid #ccc}}table{{border-collapse:collapse}}</style>
    <h1>Cas12a crRNA 设计结果</h1><p>设计工具：CaSilico casV_A；PAM：TTTV；spacer：20 nt。</p>
    <p>物种/样本：{html.escape(', '.join(s.get('name', s['id']) for s in samples))}</p>
    <p>独立组装：{len(samples)}；CaSilico 候选：{len(rows)}；导出候选：{len(selected)}。</p>
    <p>{scope}</p>
    <p>{'仅验证单一参考基因组，不能据此声称覆盖整个物种。' if len(samples)==1 else '覆盖率仅针对本次下载的组装。'}</p>
    <p>背景检查：{f'检查至多 {max_mismatches} 个碱基替换的 spacer+TTTV 匹配；不包含插入缺失、非典型 PAM 或实验交叉反应。' if background else '未提供背景数据库，特异性尚未评估。'}</p>
    <p>完整 crRNA 使用 CaSilico 默认 direct repeat；须核对所用 Cas12a 蛋白的兼容性。所有结果均为计算候选。</p>
    <table><tr>{''.join('<th>'+h+'</th>' for h in headings)}</tr>{body}</table>
    <p><a href='crRNA_designs.tsv'>完整结果</a> · <a href='crRNA_spacers.fasta'>spacer RNA</a> ·
    <a href='crRNA_full_CaSilico_default.fasta'>完整 RNA</a> · <a href='target_sites.tsv'>基因组位点</a></p>"""
    (run / "结果.html").write_text(report, encoding="utf-8")
    return {"candidates": len(rows), "selected": len(selected),
            "background_screen": f"Hamming<={max_mismatches}; TTTV; no indels" if background else "not_evaluated"}


def default_config():
    config = {"storage": str(storage_path()), "rscript": "F:/CaSilico/R/bin/Rscript.exe",
            "r_library": "F:/CaSilico/library", "wsl_distribution": "Ubuntu",
            "orthofinder": "/opt/amoebacrisp/bin/orthofinder",
            "linux_env_bin": "/opt/amoebacrisp/bin",
            "braker4": "/mnt/f/AmoebaCrisp/software/BRAKER4",
            "snakemake": "/opt/amoebacrisp/bin/snakemake",
            "direct_repeat_rna": "UAAUUUCUACUAAGUGUAGAU", "gc_min": 40, "gc_max": 65,
            "amoebadb_base": ARCHIVE}
    config.update(settings())
    config["storage"] = str(storage_path())
    return config


def annotate_braker4(samples, run, config, args):
    """Official BRAKER4 Snakemake adapter; no use of legacy braker.pl as BRAKER4."""
    needed = [s for s in samples if args.reannotate or not s.get("gff")]
    if not needed:
        return
    folder = run / "braker4"
    folder.mkdir()
    root = config["braker4"]
    if os.name == "nt":
        native_root = Path("F:/" + root[len("/mnt/f/"):]) if root.startswith("/mnt/f/") else Path(root)
    else:
        native_root = Path(root)
    template = native_root / "config.ini.example"
    if not template.exists():
        raise RuntimeError("需要 BRAKER4 注释，但尚未安装 BRAKER4；详见安装说明。已有注释时无需重新注释。")
    import configparser
    parser = configparser.ConfigParser(inline_comment_prefixes=("#",), interpolation=None)
    parser.read(template, encoding="utf-8")
    if not parser.has_section("paths"):
        parser.add_section("paths")
    parser.set("paths", "samples_file", linux_path(folder / "samples.csv"))
    parser.set("paths", "tmp_dir", linux_path(Path(config["storage"]) / "temp"))
    parser.set("paths", "busco_download_path", linux_path(Path(config["storage"]) / "cache" / "busco"))
    parser.set("paths", "compleasm_download_path", linux_path(Path(config["storage"]) / "cache" / "compleasm"))
    with (folder / "config.ini").open("w", encoding="utf-8") as handle:
        parser.write(handle)
    with (folder / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_name", "genome", "protein_fasta", "busco_lineage"])
        writer.writeheader()
        for sample in needed:
            writer.writerow({"sample_name": sample["id"], "genome": linux_path(sample["genome"]),
                             "protein_fasta": linux_path(args.protein_evidence) if args.protein_evidence else "",
                             "busco_lineage": "eukaryota_odb12"})
    command = [config["snakemake"], "--snakefile", root + "/Snakefile", "--directory", linux_path(folder),
               "--cores", str(args.threads), "--use-singularity", "--singularity-prefix",
               linux_path(Path(config["storage"]) / "software" / "containers"),
               "--singularity-args", "-B /mnt/f -B /mnt/g"]
    command = linux_command(command[0], command[1:], config)
    external(command, folder / "braker4.log", cwd=folder)
    for sample in needed:
        result_dir = folder / "output" / sample["id"] / "results"
        matches = list(result_dir.glob("braker.gff3*"))
        if len(matches) != 1:
            raise RuntimeError(f"BRAKER4 输出 GFF3 不明确：{sample['id']}；请核对当前版本的实际输出")
        sample["gff"] = str(matches[0])
        sample["annotation_method"] = "BRAKER4"


def run_pipeline(args):
    config = default_config()
    config_path = Path(args.config)
    if config_path.exists():
        config.update(json.loads(config_path.read_text(encoding="utf-8-sig")))
    storage = Path(config["storage"])
    for dirname in ["runs", "cache", "temp"]:
        (storage / dirname).mkdir(parents=True, exist_ok=True)
    os.environ.update(TEMP=str(storage / "temp"), TMP=str(storage / "temp"))
    if not Path(config["rscript"]).exists():
        raise RuntimeError("找不到 CaSilico 所需 Rscript，请检查 config.json")
    run = storage / "runs" / (dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f") + "_" + slug(args.name))
    run.mkdir()
    write_json(run / "request.json", {**vars(args), "config": config})
    log(f"输出目录：{run}")
    try:
        if args.manifest:
            samples = json.loads(Path(args.manifest).read_text(encoding="utf-8-sig"))
            for sample in samples:
                if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample["id"]):
                    raise ValueError("样本 ID 只能包含英文字母、数字、下划线、点和横线")
                for key in ["genome", "gff"]:
                    if sample.get(key):
                        sample[key] = str(Path(sample[key]).resolve())
        elif args.database == "ncbi":
            samples = ncbi_download(args.name, storage / "cache", args.max_genomes)
        else:
            samples = amoebadb_download(args.name, storage / "cache", args.max_genomes, config["amoebadb_base"])
        if not samples or len({s["id"] for s in samples}) != len(samples):
            raise ValueError("样本为空或样本 ID 重复")
        annotate_braker4(samples, run, config, args)
        write_json(run / "samples.json", samples)
        models = {}
        for sample in samples:
            models.update(extract_models(sample, run / "proteins"))
        write_json(run / "gene_map.json", {k: {a:b for a,b in v.items() if a not in ["dna", "protein"]} for k,v in models.items()})
        groups = orthogroups(samples, models, run, config, args.threads)
        groups = [(group, members) for group, members in groups
                  if all(args.min_locus <= len(models[key]["dna"]) <= args.max_locus for key in members)]
        groups.sort(key=lambda group: (max(len(models[key]["dna"]) for key in group[1]), group[0]))
        total = len(groups)
        if args.max_groups:
            groups = groups[:args.max_groups]
        write_json(run / "groups.json", {"eligible": total, "processed": len(groups), "groups": groups})
        if not groups:
            raise RuntimeError("未找到符合条件的基因组区域；不生成虚假的成功结果")
        log(f"本次设计 {len(groups)} / {total} 个可用区域；最终设计工具为 CaSilico。")
        candidates = []
        for i, (group, members) in enumerate(groups, 1):
            log(f"CaSilico：{i}/{len(groups)} {group}")
            rows = design(group, members, models, run, config, args.threads)
            candidates.extend({**row, "group": group} for row in rows)
            write_json(run / "casilico_candidates.json", candidates)
        stats = summarize(candidates, samples, run, config, args.background, args.background_mismatches)
        write_json(run / "status.json", {"status": "completed", "samples": len(samples),
                   "eligible_groups": total, "processed_groups": len(groups), **stats})
        log(f"完成：{stats['selected']} 条通过目标过滤的候选；报告：{run / '结果.html'}")
        return run
    except Exception as exc:
        write_json(run / "status.json", {"status": "failed", "error": str(exc)})
        raise


def argument_parser():
    parser = argparse.ArgumentParser(description="输入名称 → NCBI/AmoebaDB → OrthoFinder → CaSilico Cas12a crRNA")
    parser.add_argument("--name", required=True, help="拉丁学名，例如 Naegleria gruberi")
    parser.add_argument("--database", choices=["ncbi", "amoebadb"], default="ncbi")
    parser.add_argument("--config", default=str(ROOT / "config.json"))
    parser.add_argument("--manifest", help="已有本地样本 JSON，可跳过下载")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-genomes", type=int, default=0, help="0=全部独立组装")
    parser.add_argument("--max-groups", type=int, default=0, help="0=全部符合条件的区域；首轮建议填 1")
    parser.add_argument("--min-locus", type=int, default=200)
    parser.add_argument("--max-locus", type=int, default=10000)
    parser.add_argument("--background", action="append", default=[], help="背景基因组 FASTA（可重复）")
    parser.add_argument("--background-mismatches", type=int, choices=range(5), default=3,
                        help="背景允许替换的碱基数，默认3；不包含插入缺失")
    parser.add_argument("--reannotate", action="store_true", help="即使已有注释，也用 BRAKER4 重注释")
    parser.add_argument("--protein-evidence", help="BRAKER4 的多物种蛋白证据 FASTA")
    return parser


def main():
    args = argument_parser().parse_args()
    if args.threads < 1 or args.max_genomes < 0 or args.max_groups < 0 or not 24 <= args.min_locus <= args.max_locus:
        raise SystemExit("线程数、样本数或区域长度设置无效")
    run_pipeline(args)


if __name__ == "__main__":
    main()
