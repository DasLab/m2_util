# m2_util
Utilities to extra mutate-and-map 2D maps from big library sequencing runs

## TODO 

 - [ ] output to HDF5 file with [Nseq x Nres x Nres] tensor, instead of .txt file
 - [ ] provide M2 example
 - [ ] provide MOHCA example
 - [ ] convert MATLAB visualization steps into python and provide example notebook.

## Workflow

### Get the counts with `bam_to_m2.py`
First run `bam_to_m2.py`:

```
usage: bam_to_m2.py [-h] [-b [BAM ...]] -s SEQUENCES_FASTA [-mq MAP_QUALITY] [--mutdel_cutoff MUTDEL_CUTOFF] [--trim_cutoff TRIM_CUTOFF] [-o OUT_TAG]
                    [-n CHUNK_SIZE] [--start_idx START_IDX] [--end_idx END_IDX]

Processes BAM file to create 1D and 2D mut/del profiles

options:
  -h, --help            show this help message and exit
  -b, --bam [BAM ...]   BAM-formatted text file, Bowtie2
  -s, --sequences_fasta, --fasta SEQUENCES_FASTA
                        Fasta file with DNA/RNA sequences used for Bowtie2
  -mq, --map_quality MAP_QUALITY
                        minimum Bowtie2 MAPQ to consider read
  --mutdel_cutoff MUTDEL_CUTOFF
                        Filter for maximum number of mut/del in read (default 0 means no filter)
  --trim_cutoff TRIM_CUTOFF
                        Filter for maximum number of missing terminal residues (trim) in read (default 50; 0 means no filter)
  -o, --out_tag, --out_file_tag OUT_TAG
                        Tag for outfile [default: derive from first bamfile]
  -n, --chunk_size CHUNK_SIZE
                        split with this number of sequences per chunk
  --start_idx START_IDX
                        only do the reference sequences from start_idx onwards [default all]
  --end_idx END_IDX     only do the reference sequences up to end_idx [default all]

Extremely basic processing of CIGAR string in bowtie2 BAM file to produce alignments and then coding.
```

The output for each input bam will be 9 files, with suffixes: `counts_del_del.txt.gz counts_del_mut.txt.gz counts_del.txt.gz counts_mutdel_mutdel.txt.gz counts_mut_del.txt.gz counts_mutdel.txt.gz counts_mut_mut.txt.gz counts_mut.txt.gz coverage.txt.gz`

The files are:  
 - `counts_del_del.txt.gz`, ... `counts_mutdel_mutdel.txt.gz` are 2D files. These hold matrices in plain text of integer counts. The count how many reads have events -- `del` (deletion), `mut` (mutation), or `mutdel` (either deletion or mutation) -- co-occuring at residue i and j.  
	 - These 2D files take the Nres x Nres matrix for each probed reference sequence and concatenate them. The files have Nres values per row and the number of rows is Nres x Nseq, i.e., Nres rows for each reference sequence.
 - `counts_del.txt.gz`,...`counts_mut.txt.gz` are 1D files. They hold number of reads in which each event occurs at position i. The files have Nres values per row, and the number of rows is Nseq, one for each sequence.
 - `coverage.txt.gz` hold the number of reads that contributed counts to each position of each sequence. This file has Nres values per row, and the number of rows is Nseq, one for each sequence.


**Tips**:  
 - Bam file should be sorted!  
 - If you have a huge file, e.g an Ultima-generated .cram or .bam with 1B reads, you can parallelize by using `--start_idx` and `--end_idx` to select out chunks to run on different processes.  
 - If you have a bunch of `.bam` files, e.g., in different directories from a pre-split UBR run, you can supply them all with `-b` and the outputs will be in the same directories as the bam file. You can then combine the files with `merge_m2.py` (type `merge_m2.py -h` for help).
 - Merging a bunch of M2 files can get memory intensive! If not using `--start_idx`/`--end_idx`, you can provide `--chunk_size 100` to get more manageable files, and then after merging the files chunk-wise across multiple directories, concatenate them with `cat [chunk1] [chunk2]`.

 
### Visualize the maps
Once you've got the output of `bam_to_m2.py`, you can fire up MATLAB and run:

```
% Read in data
[all_counts_3d, all_counts, all_coverage, tags] = read_m2_output( filedir, data_type );

% derive maps from raw counts. 
all_cov2d = get_all_cov2d( all_counts_3d, all_counts, all_coverage [, BLANK_OUT5, BLANK_OUT3] );

% show maps with highest coverage first
coverage = get_coverage( all_coverage );
colorscale = [-5,5];
[~,sortidx] = sort(coverage,'descend');
show_m2_maps( all_cov2d, sortidx, colorscale [, BLANK_OUT5, BLANK_OUT3, coverage, headers, tags] ); 
```

**Tips**:  
- For `data_type`, use `mutdel` to capture both mutations and deletions in mutate-and-map experiments correlating SHAPE or DMS signals to intrinsic mutations and deletions or to other SHAPE/DMS events. Use `del` for MOHCA-MaP data, where it seems that only deletions carry signal.  
- `cov2d` maps are `P(i,j)/P(i)P(j) - 1`, i.e. the 2D signal above the background from 1D events.  
- An alternative is to use `get_all_m2(...)`. In that case, `colorscale = [0,0.05]` works better in `show_m2_maps`.
 
 
 
