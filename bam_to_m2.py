#! /usr/bin/env python3
import argparse
import os
import shutil
import gzip
import re
import tempfile
import shutil

parser = argparse.ArgumentParser(
                    prog = 'bam_to_2d.py',
                    description = 'Processes BAM file to create 1D and 2D mut/del profiles',
                    epilog = 'Extremely basic processing of CIGAR string in bowtie2 BAM file to produce alignments and then coding.')

parser.add_argument('-b','--bam', nargs='*',default='RTB010_CustomArray_PK50_1M7_BindTheFivePrimeEnd.bam.txt',help='BAM-formatted text file, Bowtie2')
parser.add_argument('-s','--sequences_fasta','--fasta', required=True, default='pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa',help='Fasta file with DNA/RNA sequences used for Bowtie2')
parser.add_argument('-mq','--map_quality',default=10,type=int,help='minimum Bowtie2 MAPQ to consider read')
parser.add_argument('--mutdel_cutoff',type=int,default=10,help='Filter for maximum number of mut/del in read (default 0 means no filter)' )
parser.add_argument('--trim_cutoff',type=int,default=50,help='Filter for maximum number of missing terminal residues (''trim'') in read (default 50; 0 means no filter)' )
parser.add_argument('-o','--out_tag','--out_file_tag',type=str,default='',help='Tag for outfile [default: derive from first bamfile]' )
parser.add_argument('-n','--chunk_size', default=0, type=int, help='split with this number of sequences per chunk')
parser.add_argument('--start_idx', default=0, type=int, help='only do the reference sequences from start_idx onwards [default all]')
parser.add_argument('--end_idx', default=0, type=int, help='only do the reference sequences up to end_idx [default all]')

assert(shutil.which('samtools') )
args = parser.parse_args()

def read_fasta(fasta_file):
    lines = open(fasta_file).readlines()
    sequences = []
    headers = []
    header = None
    sequence = ''
    for line in lines:
        if len(line) > 0 and line[0] == '>':
            if header is not None:
                headers.append(header)
                sequences.append(sequence)
            sequence = ''
            header = line[1:].strip('\n')
            continue
        sequence = sequence + line.strip('\n')
    if header is not None:
        headers.append(header)
        sequences.append(sequence)
    assert(len(sequences) == len(headers))
    return (sequences, headers)

(ref_sequences, ref_headers) = read_fasta(args.sequences_fasta)
ref_headers = [header.split()[0] for header in ref_headers]
ref_sequences = [sequence.replace('U','T') for sequence in ref_sequences] # need to convert to DNA!
Nref = len(ref_headers)

if len(args.bam)>1 and len(args.out_tag)>1:
    print('With multiple bam files, cannot supply --out_tag. Rerun without --out_tag and output files will have tags based on .bam files')

chunk_size = args.chunk_size
start_idx = args.start_idx
end_idx = args.end_idx
assert(not(chunk_size>0 and start_idx>0))

# Set end_idx to Nref if not specified
if end_idx == 0:
    end_idx = Nref

def initialize_counts_2d(counts_2d, coverage, Nres):
    for tag1 in ['mut','del']:
        for tag2 in ['mut','del']:
            tag = tag1 + '_' + tag2
            counts_2d[tag] = []
            for n1 in range(Nres): counts_2d[tag].append([0]*Nres)
    if len(coverage)==0:
        for n in range(Nres): coverage.append(0)
    for n in range(Nres): coverage[n] = 0

def output_to_index(ref_idx, counts_2d, coverage, fids):
    if ref_idx >= Nref or ref_idx >= end_idx:
        return (ref_idx, None)

    Nres = len(ref_sequences[ref_idx])
    fid_1d = fids[0]
    fid_2d = fids[1]

    counts_1d = {}
    for tag in ['mut','del']:
        counts_1d[tag] = [counts_2d[tag+'_'+tag][n][n] for n in range(Nres)]
        fid_out = fid_1d[tag]
        fid_out.write(",".join(map(str, counts_1d[tag])) + "\n")

    for tag1 in ['mut','del']:
        for tag2 in ['mut','del']:
            tag = tag1 + '_' + tag2
            fid_out = fid_2d[tag]
            for n1 in range(Nres):
                fid_out.write(",".join(map(str, counts_2d[tag][n1])) + "\n")

    fid_out = fid_1d['coverage']
    fid_out.write(",".join(map(str, coverage)) + "\n")

    counts_2d_all = []
    for n in range(Nres): counts_2d_all.append([0]*Nres)
    fid_out = fid_2d['mutdel']
    for n1 in range(Nres):
        for n2 in range(Nres):
            for tag1 in ['mut','del']:
                for tag2 in ['mut','del']:
                    counts_2d_all[n1][n2] += counts_2d[tag1+'_'+tag2][n1][n2]
        fid_out.write(",".join(map(str, counts_2d_all[n1])) + "\n")

    counts_1d['mutdel'] = [counts_1d['mut'][n]+counts_1d['del'][n] for n in range(Nres)]
    fid_out = fid_1d['mutdel']
    fid_out.write(",".join(map(str, counts_1d['mutdel'])) + "\n")

    ref_idx += 1
    if ref_idx < end_idx:
        initialize_counts_2d(counts_2d, coverage, len(ref_sequences[ref_idx]))

    ref_sequence = None
    if ref_idx < len(ref_sequences) and ref_idx < end_idx:
        ref_sequence = ref_sequences[ref_idx]
    return (ref_idx, ref_sequence)

def update_out_files(fids, out_tag, ref_idx, chunk_size, Nref):
    if len(fids) > 0:
        for fid_set in fids:
            for fid in fid_set.values():
                fid.close()
        fids.clear()

    num_digits = len(str(Nref))
    suffix = f'.{start_idx:0{num_digits}d}_{end_idx:0{num_digits}d}'

    fid_1d = {}
    fid_2d = {}
    for tag in ['mut', 'del', 'mutdel']:
        fid_1d[tag] = gzip.open(f'{out_tag}{suffix}.counts_{tag}.txt.gz', 'wt')

    # Separate naming for coverage file
    fid_1d['coverage'] = gzip.open(f'{out_tag}{suffix}.coverage.txt.gz', 'wt')

    for tag1 in ['mut', 'del']:
        for tag2 in ['mut', 'del']:
            tag = f'{tag1}_{tag2}'
            fid_2d[tag] = gzip.open(f'{out_tag}{suffix}.counts_{tag}.txt.gz', 'wt')

    fid_2d['mutdel'] = gzip.open(f'{out_tag}{suffix}.counts_mutdel_mutdel.txt.gz', 'wt')

    fids.append(fid_1d)
    fids.append(fid_2d)

# Create regions file for samtools
regions_tag = ''
if start_idx > 0:
    temp_dir = tempfile.mkdtemp()
    regions_file = f'{temp_dir}/{args.sequences_fasta}.{start_idx:07d}_{end_idx:07d}.bed'
    with open(regions_file, 'w') as fid_regions:
        for i in range(start_idx-1, end_idx):
            fid_regions.write(f'{ref_headers[i]} 0 {len(ref_sequences[i])}\n')
    regions_tag = f' --regions-file {regions_file}'

for bam in args.bam:
    print(f'Doing BAM file {bam}')
    assert(bam.endswith('.bam') or bam.endswith('.cram'))

    out_tag = args.out_tag
    if len(out_tag)==0: out_tag = bam.replace('.bam','').replace('.cram','')

    command = f"samtools view {bam} {regions_tag} | awk '$5>={args.map_quality}'"
    fid_bam = os.popen(command)

    line = fid_bam.readline()
    cols = line.split()
    ref_idx = start_idx - 1
    header = ''

    count = 0
    count_filter = 0

    fids = []

    header = line.split()[2]
    ref_sequence = ref_sequences[ref_headers.index(header)]
    Nres = len(ref_sequence)
    update_out_files(fids, out_tag, ref_idx, chunk_size, Nref)
    counts_2d = {}
    coverage = []
    initialize_counts_2d(counts_2d, coverage, Nres)

    while line:
        cols = line.split()
        next_header = cols[2]
        if header != next_header:
            while ref_headers[ref_idx] != next_header:
                ref_idx, ref_sequence = output_to_index(ref_idx, counts_2d, coverage, fids)
                if ref_sequence is None:
                    break
            header = next_header
            if ref_sequence is None:
                break

        start_pos = int(cols[3])
        cigar = cols[5]
        signed_tmpl_len = int(cols[8])
        read = cols[9]

        seqa = ''
        read_pos = 0
        seqa += '.' * (start_pos - 1)

        # Parse CIGAR with regex: [('22', 'M'), ('1', 'I'), ...]
        cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

        for length_str, op in cigar_tuples:
            length = int(length_str)
            if op == 'M':  # match or mismatch
                seqa += read[read_pos:read_pos + length]
                read_pos += length
            elif op == 'I':  # insertion to the reference — skip, not in alignment
                read_pos += length
            elif op == 'D':  # deletion from the reference — gap in read
                seqa += '-' * length
            elif op == 'S':  # soft clip — skip bases in read
                read_pos += length
            elif op == 'H':  # hard clip — ignore entirely
                continue
            elif op == 'N':  # skipped region — typically introns, same as D here
                seqa += '-' * length
            # P (padding) and others are rarely used; can be ignored or handled as needed

        assert(len(seqa) <= len(ref_sequence))
        end_pos = len(seqa)

        if signed_tmpl_len < 0:
            seqa = '.' * (len(ref_sequence) - len(seqa)) + seqa
        else:
            seqa = seqa + '.' * (len(ref_sequence) - len(seqa))

        assert len(seqa) == len(ref_sequence)

        pos = {'mut': [], 'del': []}
        for n, (ref_nt, read_nt) in enumerate(zip(ref_sequence, seqa)):
            if ref_nt not in 'ACGT-': continue
            if read_nt not in 'ACGT-': continue
            if read_nt == '-':  pos['del'].append(n)
            elif ref_nt != read_nt: pos['mut'].append(n)

        count += 1
        if count % 1000 == 0: print(f'Processed {count:10d} lines...  and number that pass filter: {count_filter:10d}')

        total_mutdel = len(pos['mut']) + len(pos['del'])
        total_trim = seqa.count('.')
        if (args.mutdel_cutoff == 0 or total_mutdel <= args.mutdel_cutoff) and \
           (args.trim_cutoff == 0 or total_trim <= args.trim_cutoff):
            count_filter += 1
            for tag1 in ['mut','del']:
                for tag2 in ['mut','del']:
                    for n1 in pos[tag1]:
                        for n2 in pos[tag2]:
                            counts_2d[tag1+'_'+tag2][n1][n2] += 1
            for n in range(start_pos-1, end_pos):
                coverage[n] += 1

        line = fid_bam.readline()

    while ref_idx < end_idx:
        ref_idx, ref_sequence = output_to_index(ref_idx, counts_2d, coverage, fids)

    for fid_set in fids:
        for fid in fid_set.values():
            fid.close()

    print(f'Processed {count} reads. Filtered {count_filter} reads with <= {args.mutdel_cutoff} mutations+deletions and <= {args.trim_cutoff} trimmed nts at termini.')

if len(regions_tag)>0: shutil.rmtree(temp_dir)



