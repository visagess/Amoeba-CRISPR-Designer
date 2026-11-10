"""Extract explicitly annotated nuclear 18S loci from the public archive."""
import gzip
import hashlib
from pathlib import Path
import re
from urllib.parse import urlencode
import xml.etree.ElementTree as ET
import amoebacrisp as app
from ssu18s_import import read_records


def retrieve_features(root,target_taxon,genera,min_bp,max_bp):
    from ssu18s import EUTILS,search_ids,identifies_18s,classify
    rows=[];rejected=[];queries={};accessions=set()
    config=app.default_config()
    for genus in genera:
        ids=search_ids('taxonomy',f'txid{genus["taxid"]}[Subtree] AND species[Rank]')
        species=[]
        for start in range(0,len(ids),100):
            xml=ET.fromstring(app.fetch(EUTILS+'efetch.fcgi?'+urlencode({'db':'taxonomy','id':','.join(ids[start:start+100])})))
            species.extend(t.findtext('ScientificName') for t in xml.findall('Taxon') if t.findtext('ScientificName'))
        samples=app.amoebadb_download(genus['name'],root/'archive_cache',0,config['amoebadb_base'])
        queries[genus['name']]=[s['id'] for s in samples]
        for sample in samples:
            accessions.add(sample['id'])
            code=sample['id'].split('_',1)[1]
            matches=[s for s in species if len(s.split())==2 and code.casefold().startswith((s[0]+s.split()[1]).casefold())]
            if matches:
                longest=max(len(s.split()[1]) for s in matches)
                matches=[s for s in matches if len(s.split()[1])==longest]
            if len(matches)!=1 or not sample.get('gff'):
                rejected.append({'accession':sample['id'],'organism':'','reason':'archive_species_unresolved_or_no_annotation'})
                continue
            organism=matches[0]
            genome={header.split()[0]:(sequence,header) for header,sequence in read_records(sample['genome'])}
            opener=gzip.open if str(sample['gff']).endswith('.gz') else open
            features=[];organelle_contigs=set()
            with opener(sample['gff'],'rt',encoding='utf-8') as handle:
                for line in handle:
                    if line.startswith('##FASTA'):break
                    if line.startswith('#'):continue
                    fields=line.rstrip().split('\t')
                    if len(fields)!=9:continue
                    attrs=app.gff_attributes(fields[8])
                    label='; '.join(attrs.get(k,'') for k in ('product','gene','Name','standard_name','description'))
                    if fields[2].casefold() in ('region','chromosome','contig') and re.search('mitochond|plastid|chloroplast',fields[8],re.I):
                        organelle_contigs.add(fields[0])
                    if fields[2].casefold() not in ('rrna','misc_rna','gene'):continue
                    if identifies_18s(label) and not re.search(r'internal transcribed|5[.]8s|large subunit|28s',label,re.I):
                        features.append((fields,attrs,label))
            grouped={}
            for fields,attrs,label in features:
                key=attrs.get('ID') or ':'.join(fields[:1]+fields[3:7])
                grouped.setdefault(key,[]).append((fields,attrs,label))
            loci=[]
            for parts in grouped.values():
                if len({(f[0],f[6]) for f,a,l in parts})!=1:
                    rejected.append({'accession':sample['id'],'organism':organism,'reason':'multi_contig_or_strand_18S_annotation'})
                    continue
                fields,attrs,label=parts[0];fields=fields.copy()
                fields[3]=str(min(int(f[3]) for f,a,l in parts));fields[4]=str(max(int(f[4]) for f,a,l in parts))
                loci.append((fields,attrs,label))
            seen=set()
            for fields,attrs,label in loci:
                contig=fields[0];begin,end=int(fields[3]),int(fields[4]);strand=fields[6]
                key=(contig,begin,end,strand)
                if key in seen:continue
                seen.add(key)
                if contig not in genome or not 1<=begin<=end<=len(genome[contig][0]) or strand not in ('+','-'):
                    rejected.append({'accession':sample['id'],'organism':organism,'reason':'invalid_18S_coordinates'})
                    continue
                sequence,title=genome[contig]
                if contig in organelle_contigs or re.search('mitochond|plastid|chloroplast',title+' '+fields[8],re.I):continue
                dna=sequence[begin-1:end]
                if strand=='-':dna=app.rc(dna)
                accession=sample['id']+':'+contig
                rows.append({'id':app.slug(sample['id']+'_'+contig)+f'__18S_{begin}_{end}_{strand}',
                             'accession':accession,'organism':organism,'role':classify(organism,target_taxon['name'],[g['name'] for g in genera]),
                             'title':label,'taxid':'','strain':'','isolate':'','sequence':dna,'length':len(dna),
                             'source_url':sample['url'],'source_record_length':len(sequence),'DNA_supported':True,
                             'mol_type':'genomic DNA','annotated_partial':attrs.get('partial','false').casefold()=='true',
                             'ambiguous_bases':sum(b not in 'ACGT' for b in dna),'sha256':hashlib.sha256(dna.encode()).hexdigest(),
                             'record_start_1based':begin,'record_end_1based':end,'strand':strand,
                             'genomic_introns_preserved':True,'annotation_method':'AmoebaDB_GFF3_explicit_18S',
                             'archive_release':sample['release']})
    return rows,rejected,queries,accessions
