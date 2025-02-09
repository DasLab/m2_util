#! /usr/bin/env python3
import argparse
import os
import shutil
import gzip

parser = argparse.ArgumentParser(
                    prog = 'bam_to_2d.py',
                    description = 'Processes BAM file to create 1D and 2D mut/del profiles',
                    epilog = 'Extremely basic processing of CIGAR string in bowtie2 BAM file to produce alignments and then coding.')

parser.add_argument('-b','--bam', nargs='*',default='RTB010_CustomArray_PK50_1M7_BindTheFivePrimeEnd.bam.txt',help='BAM-formatted text file, Bowtie2')
parser.add_argument('-s','--sequences_fasta','--fasta', required=True, default='pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa',help='Fasta file with DNA/RNA sequences used for Bowtie2')
parser.add_argument('-mq','--map_quality',default=10,type=int,help='minimum Bowtie2 MAPQ to consider read')
parser.add_argument('--mut_cutoff',type=int,default=0,help='Filter for maximum number of mut/del in read (default 0 means no filter)' )
parser.add_argument('-o','--out_tag','--out_file_tag',type=str,default='',help='Tag for outfile [default: derive from first bamfile]' )
parser.add_argument('-n','--chunk_size', default=0, type=int, help='split with this number of sequences per chunk')

assert(shutil.which('samtools') )
args = parser.parse_args()

def read_fasta( fasta_file ):
    lines = open( fasta_file ).readlines()
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
    assert( len(sequences) == len(headers ) )
    return (sequences,headers)

(ref_sequences,ref_headers) = read_fasta( args.sequences_fasta )
ref_headers = [header.split()[0] for header in ref_headers]
Nref = len(ref_headers)

out_tag = args.out_tag
if len( out_tag ) == 0:
    if len( args.bam ) == 1:
        out_tag = args.bam[0].replace('.bam.txt','').replace('.bam','').replace('.bam','')
    else:
        print('Supply out_tag with -o since you have multiple bam files.')
        exit(0)

def initialize_counts_2d( counts_2d, coverage, Nres ):
    for i in ['mut','del']:
        if i not in counts_2d: counts_2d[i] = {}
        for j in ['mut','del']:
            counts_2d[i][j] = []
            for n in range(Nres):
                counts_2d[i][j].append( [0]*Nres )
    if len(coverage)==0:
        for n in range(Nres): coverage.append(0)
    for n in range(Nres): coverage[n] = 0

def update_out_files( fids, out_tag,ref_idx, chunk_size, Nref, gzipped=True ):
    if gzipped:
        file_open = gzip.open
        gzip_ext=".gz"
    else:
        file_open=open
        gzip_ext=""
    if len(fids)>0:
        for (fid_2d,fid_1d,fid_coverage) in fids[-1]:
            for tag1 in ['mut','del']:
                fid_1d[tag1].close()
                for tag2 in ['mut','del']:
                    fid_2d[tag1][tag2].close()
            fid_2d['mutdel'].close()
            fid_1d['mutdel'].close()
            fid_coverage.close()

    split_tag=''
    if chunk_size>0 and Nref>chunk_size:
        chunk_start = ref_idx+1
        chunk_end   = min( ref_idx + chunk_size, len(ref_headers))
        split_tag='.%07d_%07d' % (chunk_start,chunk_end)
    if ref_idx < Nref-1:
        fid_2d = {}
        fid_1d = {}
        for tag1 in ['mut','del']:
            fid_2d[tag1] = {}
            for tag2 in ['mut','del']:
                outfile = out_tag + '.counts_%s_%s.txt' % (tag1,tag2)
                fid_2d[tag1][tag2] = file_open( outfile+gzip_ext, 'wt')
            outfile = out_tag + '.counts_%s.txt' % (tag1)
            fid_1d[tag1] = file_open( outfile+gzip_ext, 'wt' )
        outfile = out_tag + '.counts_%s_%s.txt' % ('mutdel','mutdel')
        fid_2d['mutdel'] = file_open( outfile+gzip_ext, 'wt' )
        outfile = out_tag + '.counts_%s.txt' % ('mutdel')
        fid_1d['mutdel'] = file_open( outfile+gzip_ext, 'wt' )
        outfile = out_tag + '.coverage.txt'
        fid_1d['coverage'] = file_open( outfile+gzip_ext, 'wt' )
        fids.append( (fid_2d,fid_1d) )

def output_to_index( ref_idx, counts_2d, coverage, fids ):
    (fid_2d, fid_1d) = fids[-1]
    Nres = len(counts_2d['mut']['mut'])
    # 2D profiles
    for tag1 in ['mut','del']:
        for tag2 in ['mut','del']:
            for n1 in range(Nres):
                fid_2d[tag1][tag2].write(",".join(map(str, counts_2d[tag1][tag2][n1])) + "\n")

    # 1D profiles
    counts_1d={}
    for tag in ['mut','del']:
        counts_1d[tag] = [counts_2d[tag][tag][n][n] for n in range(Nres)]
        fid_1d[tag].write( ",".join(map(str,counts_1d[tag]))+"\n")
    fid_1d['coverage'].write( ",".join(map(str,coverage))+"\n")

    # combine mut/del
    counts_2d_all = []
    for n in range( Nres ): counts_2d_all.append( [0]*Nres )
    fid_out = fid_2d['mutdel']
    for n1 in range(Nres):
        for n2 in range(Nres):
            for tag1 in ['mut','del']:
                for tag2 in ['mut','del']:
                    counts_2d_all[n1][n2] += counts_2d[tag1][tag2][n1][n2]
        fid_out.write(",".join(map(str, counts_2d_all[n1])) + "\n")

    counts_1d['mutdel'] = [ counts_1d['mut'][n]+counts_1d['del'][n] for n in range(Nres)]
    fid_out = fid_1d['mutdel']
    fid_out.write(",".join(map(str, counts_1d['mutdel'])) + "\n")

    ref_idx += 1
    initialize_counts_2d(counts_2d, coverage, Nres)
    if chunk_size>0 and (ref_idx+1) % chunk_size == 1:
        update_out_files( fids, out_tag,ref_idx, chunk_size, Nref)

    ref_sequence = None
    if ref_idx < len(ref_sequences): ref_sequence = ref_sequences[ref_idx]
    return (ref_idx, ref_sequence)

chunk_size = args.chunk_size
for bam in args.bam:

    print( 'Doing BAM file %s' % bam )
    assert( bam[-4:]=='.bam' or  bam[-4:]=='.bam' )

    out_tag = args.out_tag
    if len(out_tag)==0: out_tag = bam.replace('.bam','')

    command = "samtools view %s | awk '$5>=%d'" % (bam,args.map_quality )
    fid_bam = os.popen( command )

    line = fid_bam.readline()
    cols = line.split()
    ref_idx = 0 # where we are in reference sequence file, 0 indexed
    header = ''

    count = 0
    count_filter = 0

    fids = [] # outfiles

    header = line.split()[2]
    ref_sequence = ref_sequences[ ref_headers.index( header ) ]
    Nres = len(ref_sequence)
    update_out_files( fids, out_tag,ref_idx, chunk_size, Nref)
    counts_2d = {}
    coverage = []
    initialize_counts_2d( counts_2d, coverage, Nres )

    while line:
        cols = line.split()
        next_header = cols[2]
        if header != next_header: # we're ready to output the 2D stats
            while ref_headers[ref_idx] != next_header:
                ref_idx, ref_sequence = output_to_index( ref_idx, counts_2d, coverage, fids )
        start_pos = int(cols[3])
        cigar = cols[5]
        signed_tmpl_len = int(cols[8])
        read = cols[9]

        #########################################
        # create alignment output seqa for read.
        #########################################
        # parse cigar string, e.g., "22M1I23M1D69M" => "22M 1I 23M 1D 69M"
        cpos = 0 # cigar position
        spos = 0 # sequence position
        seqa = ''
        seqa += '.'*(start_pos-1)
        for k,s in enumerate(cigar):
            num = cigar[cpos:(k+1)]
            if not num.isnumeric():
                nres = int(cigar[cpos:k])
                indelcode = cigar[k]
                if indelcode == 'M':
                    seqa += read[spos:(spos+nres)]
                    spos += nres
                elif indelcode == 'D':
                    seqa += '-'*nres
                    spos += 0
                elif indelcode == 'I':
                    spos += nres
                cpos = k+1 # advance to next chunk of cigar string
        assert(len(seqa) <= len(ref_sequence))
        end_pos = len(seqa)
        if signed_tmpl_len < 0:
            seqa = '.'*(len(ref_sequence)-len(seqa)) + seqa
        else:
            seqa = seqa + '.'*(len(ref_sequence)-len(seqa))
        assert(len(seqa)==len(ref_sequence))

        pos = {}
        pos['mut'] = []
        pos['del'] = []
        for n,nt_pair in enumerate(zip(ref_sequence,seqa)):
            if nt_pair[0] not in 'ACGT-': continue
            if nt_pair[1] not in 'ACGT-': continue
            if nt_pair[1] == '-': pos['del'].append(n)
            elif nt_pair[0] != nt_pair[1]: pos['mut'].append(n)

        count += 1
        if count % 1000 == 0: print('Processed %d lines...' % count)

        total_mutdel = len(pos['mut']) + len(pos['del'])
        if args.mut_cutoff == 0 or total_mutdel <= args.mut_cutoff:
            count_filter += 1
            for tag1 in ['mut','del']:
                for tag2 in ['mut','del']:
                    for n1 in pos[tag1]:
                        for n2 in pos[tag2]:
                            counts_2d[tag1][tag2][n1][n2] += 1
            for n in range(start_pos-1,end_pos): coverage[n] += 1

        line = fid_bam.readline()

    while ref_idx<Nref:
        ref_idx, ref_sequence = output_to_index( ref_idx, counts_2d, coverage, fids )

    fid_bam.close()
    print('Processed %d from %d lines\n' % (count_filter,count) )
    print('Created %d sets of files with names like:' % len(fids) )
    (fid_2d, fid_1d) = fids[0]
    for tag1 in ['mut','del']:
        for tag2 in ['mut','del']:
            print(fid_2d[tag1][tag2].name)
    for tag in ['mut','del','mutdel','coverage']: print( fid_1d[tag].name)
