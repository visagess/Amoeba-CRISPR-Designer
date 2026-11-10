"""Exhaustive substitution-only Cas12a background screening using NumPy.

Every unambiguous 20 nt site immediately following a canonical TTTV PAM is
encoded on both strands. Pairwise two-bit XOR/popcount gives the exact Hamming
distance. This has no seed or local-alignment sensitivity gap.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import numpy as np

from amoebacrisp import checksum, read_fasta, rc, write_json

DTYPE = np.dtype([('code','<u8'), ('contig','<u4'), ('start','<u4'), ('strand','i1'), ('pam','u1')])
LOOKUP = np.full(256, 255, dtype=np.uint8)
for index, letter in enumerate(b'ACGT'):
    LOOKUP[letter] = index
PAIR_MASK = np.uint64(0x5555555555)
BYTE_POPCOUNT = np.array([int(i).bit_count() for i in range(256)], dtype=np.uint8)


def popcount(values):
    if hasattr(np, 'bitwise_count'):
        return np.bitwise_count(values)
    raw = np.ascontiguousarray(values).view(np.uint8).reshape(-1,8)
    return BYTE_POPCOUNT[raw].sum(axis=1)


def code_for(sequence):
    if not re.fullmatch('[ACGT]{20}',sequence):
        raise ValueError('需要 20 nt 的明确 DNA spacer')
    result=0
    for base in sequence:
        result=(result<<2) | 'ACGT'.index(base)
    return np.uint64(result)


def encode_sites(sequence, contig):
    if len(sequence) >= 2**32:
        raise ValueError('单条 contig 超出本版 32-bit 坐标范围')
    for strand, dna in [(1,sequence),(-1,rc(sequence))]:
        if len(dna)<24:
            continue
        values=LOOKUP[np.frombuffer(dna.encode('ascii'),dtype=np.uint8)]
        positions=np.flatnonzero((values[:-23]==3) & (values[1:-22]==3) &
                                 (values[2:-21]==3) & (values[3:-20]<3))
        codes=np.zeros(len(positions),dtype=np.uint64)
        valid=np.ones(len(positions),dtype=bool)
        for offset in range(20):
            bases=values[positions+4+offset]
            valid &= bases<4
            codes=(codes<<np.uint64(2)) | bases.astype(np.uint64)
        positions,codes=positions[valid],codes[valid]
        rows=np.empty(len(positions),dtype=DTYPE)
        rows['code']=codes
        rows['contig']=contig
        rows['start']=positions+5 if strand==1 else len(dna)-positions-23
        rows['strand']=strand
        rows['pam']=values[positions+3]
        yield rows


def build_index(background, path):
    path=Path(path)
    stamp=Path(str(path)+'.json')
    source_hash=checksum(background)
    if path.exists() and stamp.exists():
        metadata=json.loads(stamp.read_text())
        if metadata.get('background_sha256')==source_hash and path.stat().st_size==metadata.get('records',-1)*DTYPE.itemsize:
            return metadata
    temporary=Path(str(path)+'.part')
    path.parent.mkdir(parents=True,exist_ok=True)
    records,contigs=0,[]
    with temporary.open('wb') as handle:
        for number,(name,seq) in enumerate(read_fasta(background)):
            contigs.append({'id':number,'name':name,'length':len(seq)})
            for rows in encode_sites(seq,number):
                rows.tofile(handle)
                records += len(rows)
            if (number+1)%5000==0:
                print(f'PAM index: {number+1} contigs, {records} valid sites',flush=True)
    temporary.replace(path)
    metadata={'background_sha256':source_hash,'records':records,'contigs':contigs,
              'spacer_length':20,'pam':'TTTV','strands':'both','ambiguous_sites':'excluded',
              'dtype':DTYPE.descr}
    write_json(stamp,metadata)
    print(f'PAM index completed: {len(contigs)} contigs, {records} sites',flush=True)
    return metadata


def scan_queries(queries,index,metadata,output,summary,max_mismatches=3,chunk_size=2_000_000):
    queries=list(read_fasta(queries))
    packed=[(name,sequence,code_for(sequence)) for name,sequence in queries]
    count={name:{'sequence':seq,'exact_hits':0,'near_hits':0} for name,seq,_ in packed}
    total=metadata['records']
    data=np.memmap(index,dtype=DTYPE,mode='r',shape=(total,)) if total else np.empty(0,dtype=DTYPE)
    with Path(output).open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['query','spacer_dna','subject','start','end','strand','pam','mismatches'],delimiter='\t')
        writer.writeheader()
        for offset in range(0,total,chunk_size):
            chunk=data[offset:offset+chunk_size]
            values=np.ascontiguousarray(chunk['code'])
            for name,sequence,query_code in packed:
                difference=values ^ query_code
                distance=popcount((difference | (difference>>np.uint64(1))) & PAIR_MASK)
                matching=np.flatnonzero(distance<=max_mismatches)
                count[name]['near_hits'] += len(matching)
                count[name]['exact_hits'] += int(np.count_nonzero(distance==0))
                for at in matching:
                    hit=chunk[at]
                    start=int(hit['start'])
                    writer.writerow({'query':name,'spacer_dna':sequence,
                                     'subject':metadata['contigs'][int(hit['contig'])]['name'],
                                     'start':start,'end':start+19,'strand':'+' if hit['strand']==1 else '-',
                                     'pam':'TTT'+'ACG'[int(hit['pam'])],'mismatches':int(distance[at])})
            print(f'PAM screen: {min(offset+chunk_size,total)}/{total} sites, {len(packed)} spacers',flush=True)
    result={'method':'exhaustive_2bit_Hamming','background_sha256':metadata['background_sha256'],
            'max_substitutions':max_mismatches,'insertions_deletions':False,
            'pam':'TTTV','spacer_length':20,'strands':'both','ambiguous_sites':'excluded','site_count':total,'queries':count}
    write_json(summary,result)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--background',required=True)
    parser.add_argument('--queries')
    parser.add_argument('--index',required=True)
    parser.add_argument('--output')
    parser.add_argument('--summary')
    parser.add_argument('--build-only',action='store_true')
    parser.add_argument('--mismatches',type=int,default=3,choices=range(5))
    args=parser.parse_args()
    metadata=build_index(args.background,args.index)
    if args.build_only:
        return
    if not all([args.queries,args.output,args.summary]):
        parser.error('筛查时必须提供 queries、output 和 summary')
    scan_queries(args.queries,args.index,metadata,args.output,args.summary,args.mismatches)


if __name__=='__main__':
    main()
