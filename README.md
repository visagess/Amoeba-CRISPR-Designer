# Amoeba CRISPR Designer

阿米巴 CRISPR 设计。A local workflow for comparing amoeba sequences and designing candidate CRISPR-Cas12a crRNAs using CaSilico and EasyDesign.

The current primary workflow compares nuclear 18S rDNA. It finds conserved and conditionally species-specific windows, separates windows with and without a canonical TTTV PAM, and exports Excel results with GC content, species coverage, exception species and matching evidence.

## 中文快速说明

当前主流程为核18S分析，可在线获取已验证的 N. gruberi 专题数据，也可导入自己的FASTA并选择目标文件。输入长度限定800–4000 bp，使用MAFFT比对，分别输出有PAM和无PAM的保守、特异候选。特异区逐段统计例外物种，默认允许3种，未定种命中另列待核查。CaSilico保留20nt独立设计，EasyDesign对有PAM的21nt候选独立设计和预测，最后导出Excel。

首次使用请先完成下方依赖安装，将 `config.example.json` 复制为 `config.json` 并填写本机工具路径，再运行 `gui.py` 或 `launch.cmd`。本机配置文件由Git忽略。README下方给出命令行、本地导入格式、筛选规则和结果说明。

## Workflow

```text
NCBI nuclear 18S retrieval / local per-species FASTA files
  → annotation-aware extraction and 800–4000 bp length filter
  → MAFFT alignment and physical DNA coordinate projection
  → conserved / specific windows, with / without TTTV PAM
  → CaSilico 20 nt designs and background screening
  → independent EasyDesign 21 nt design and activity prediction
  → Excel workbook and underlying audit files
```

The online 18S preset is currently validated for target **Naegleria gruberi** and background genera **Naegleria + Entamoeba** through NCBI Nucleotide. The local FASTA workflow supports a user-selected target file and background files with confirmed species labels. The GUI also retains genome-region and shared-gene modes. Those modes support NCBI Datasets or the public AmoebaDB archive. BRAKER4 and OrthoFinder belong to the shared-gene mode.

## Features

- Search multiple 18S / SSU aliases and extract annotated nuclear small-subunit sequences.
- Preserve sequence provenance, different within-species variants, genomic introns and continuous DNA coordinates.
- Import multiple FASTA files, confirm one species per file, and select the target file.
- Adjust conserved-species coverage, mismatch tolerance and per-window exception-species allowance.
- Screen both strands of all supplied background sequences using a seeded exhaustive substitution search.
- Count a matching species once even when multiple records match. Report unresolved labels separately.
- Produce four comparison groups plus independently evaluated EasyDesign Cas12a designs in Excel.

## Installation

The tested deployment uses **Windows with Python 3.11 and an Ubuntu WSL environment**. New software, datasets and caches can be kept on F\ drive. Linux tools run inside WSL. To keep their physical storage on F\ drive, place the WSL distribution's virtual disk there as well.

1. Download or clone this repository into `F:\AmoebaCrisp\app`.
2. Create a Windows Python environment and install the public export dependencies.

```powershell
python -m venv F:\AmoebaCrisp\software\app-env
F:\AmoebaCrisp\software\app-env\Scripts\python.exe -m pip install --cache-dir F:\AmoebaCrisp\software\pip-cache -r requirements.txt
Copy-Item config.example.json config.json
```

3. In Ubuntu WSL, use Conda/Mamba to install the alignment and comparison tools. Run from the repository directory, adjusting the prefix to your WSL installation.

```bash
conda env create -p "$HOME/amoebacrisp-env" -f environment.yml
```

Set `linux_env_bin` in `config.json` to that environment's absolute Linux `bin` path. It must contain `python`, `mafft`, `blastn`, `makeblastdb` and `dustmasker`. Set `wsl_distribution` to the distribution shown by `wsl -l -q`.

4. Install **CaSilico** using its upstream instructions. Set `rscript` and `r_library` to your existing installation. The R library must load `CaSilico`. This project calls its local R interface.
5. Install **EasyDesign** and its official Cas12a classification and regression models from the upstream repository. Keep it in a separate Python 3.8 / TensorFlow 2.3 environment matching the tested model runtime. Set `easydesign_source` and `easydesign_python` in `config.json`. The source directory must contain `bin/design.py` and `easyDesign/models/{classify,regress}/cas12a`.

The main `requirements.txt` covers the Windows controller and Excel exporter. Third-party bioinformatics tools, R packages, model weights and the legacy EasyDesign runtime require their own upstream installation. See [DEPENDENCIES.md](DEPENDENCIES.md).

`config.json` is local configuration and is ignored by Git. `AMOEBA_CRISPR_CONFIG` can select an alternative config file. `AMOEBA_CRISPR_HOME` can override the storage root. Defaults use `F:/AmoebaCrisp`.

## Usage

### Graphical interface

```powershell
F:\AmoebaCrisp\software\app-env\Scripts\python.exe gui.py
```

Alternatively, set `AMOEBA_PYTHON` to the environment's Python executable and double-click `launch.cmd`.

For the online pilot, select `18S专题（N. gruberi）` and keep its verified target/background preset. For local data, click **导入FASTA并选目标文件**, select unaligned nucleotide FASTA files, confirm each species label, and select exactly one target file. A file may contain multiple variants of the same species. Confirm that the supplied sequences are nuclear 18S.

### Command line

```powershell
python ssu18s_run.py --threads 4 --mismatches 2 --background-min-mismatches 2 --max-exception-species 3
```

Add `--reuse-data` to reuse the curated inputs. Optional `--coverage-min` and `--coverage-max` select the conserved-species interval, including the target species. Leaving both unset requires coverage of the available same-genus species set.

A local import uses `--local-manifest path/to/inputs.json` with this structure.

```json
{
  "files": [
    {"path": "F:/inputs/target.fasta", "species": "Naegleria gruberi", "target": true},
    {"path": "F:/inputs/background.fasta", "species": "Naegleria pagei", "target": false}
  ]
}
```

The GUI constructs this manifest automatically. Each local run gets a separate data directory with copies of the source files.

### Genome modes

```powershell
python region_pipeline.py --target "Naegleria gruberi" --background-genus Naegleria --background-genus Entamoeba --database ncbi --design-regions 10
python amoebacrisp.py --name "Naegleria gruberi" --database ncbi --max-groups 10
```

Shared-gene mode uses existing annotations when available. For reannotation it calls the official BRAKER4 Snakemake workflow, which requires its own installation, containers and evidence configuration. Multi-assembly shared-gene identification requires OrthoFinder. The default 18S workflow compares the marker directly.

## Default screening rules

| Setting | Default |
|---|---|
| Extracted nuclear 18S length | 800–4000 bp, inclusive |
| Comparison window | 20 nt |
| Target / conserved tolerance | At most 2 substitutions |
| Conserved coverage | All available same-genus species, adjustable |
| Background distance | At least 2 substitutions for non-exception species |
| Allowed non-specific named species | At most 3 per window, adjustable |
| PAM | TTTV on either reference strand |
| EasyDesign spacer | 21 nt, independently screened |

All background records remain in screening. A 0–1 mismatch hit in any record makes that named species an exception under the default background threshold. Each window has its own exception list. The default cap of 3 was selected for the pilot after cap 2 produced no candidates.

`sp.`, `cf.`, `aff.` and other unresolved labels are screened separately. At cap 0, near matches to these records reject the candidate. At a positive cap, unresolved matches are reported and the candidate is conditional pending review. Missing homologous-site coverage is shown as unassessed. The distinguishable-species count includes only named species with site coverage and no near match.

Different within-species sequences remain in comparison. Identical sequences within a species are deduplicated with provenance preserved. Conserved coverage requires all covered variants of a species to fit the tolerance. Specificity screening uses any background-record match as evidence.

## Outputs

Excel workbooks are written under `F:/AmoebaCrisp/outputs` by default. The six sheets are overview, conserved with PAM, conserved without PAM, specific with PAM, specific without PAM, and EasyDesign 21 nt designs. Empty groups have headers only.

Rows include coordinates, strand, GC content, coverage, exception species and matching evidence where applicable. Underlying JSON, TSV, FASTA, alignments, CaSilico results, official EasyDesign design output and model/input checksums remain under the run's `comparison` directory.

EasyDesign scores are official model outputs for the reference target. They are not efficiency percentages. The inactive classification gives 0, otherwise the shifted activity is `max(0, regression + 4)`. The official high-activity threshold is retained. Scores from 21 nt designs are kept separate from CaSilico 20 nt designs. Complete RNA strings use the configured default direct repeat plus spacer. The PAM remains part of the DNA target.

## Validation and scope

Run the automated tests with `python -m unittest discover -s tests -v`. Tests cover extraction, local import, strand handling, continuous-coordinate projection, mismatch screening, per-window exception counting and model-context construction.

The original local deployment was exercised with real N. gruberi inputs, MAFFT, CaSilico and the official EasyDesign Cas12a model. BRAKER4 reannotation and multi-assembly OrthoFinder require additional end-to-end validation. The AmoebaDB option uses a public CyVerse archive whose versions may differ from the current database website. The online 18S preset uses NCBI Nucleotide.

The distributable was also tested end-to-end on 486 eligible input sequences. That fresh alignment produced 14 conserved PAM windows, 310 conserved no-PAM windows, 1 conditional specific PAM window, 7 conditional specific no-PAM windows, and 15 EasyDesign designs. Multithreaded MAFFT iterative refinement can vary between fresh runs, changing some conserved windows. Preserve `comparison/MAFFT/aligned.fasta` and its checksum to reproduce a particular result. The earlier cached alignment had 15/308 conserved windows and 16 EasyDesign designs. Both runs retained the same count of 8 conditional specific 20 nt windows.

Sequence-screening conclusions apply to the selected inputs and mismatch rules. Candidate windows may overlap. The 18S workflow evaluates background 18S similarity and does not establish whole-genome off-target absence. Substitution screening does not model indels or experimental mismatch tolerance. No-PAM windows remain comparison results requiring an appropriate Cas12a targeting strategy before guide use.

## Acknowledgements

The workflow requirements were inspired by [CrispHunter](https://github.com/2021JohnSheng/CrispHunter). This repository contains an independent implementation for the amoeba use case. Third-party software remains under its upstream licenses. See [DEPENDENCIES.md](DEPENDENCIES.md).
