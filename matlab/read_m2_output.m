function [all_counts_3d, all_counts, all_coverage, tags] = read_m2_output( filedir, data_type, tags );
% [all_counts_3d, all_counts, all_coverage, tags] = read_m2_output( filedir, data_type, tags );
%
% Read in counts files output by bam_to_m2. Currently looks for events of
%   one data_type, no correlations of one kind of event and another.
%
% Inputs:
%
%   filedir   = [string] directory with files like *.coverage.txt.gz,
%                 *.counts_mutdel_mutdel.txt.gz, output by bam_to_m2.py
%   data_type = [string] which data_type to look for correlations.
%                   Should be 'mutdel' (mutations and deletions in 2A3 or DMS M2) 
%                    'del' (best for MOHCA-MaP-seq), or 'mut' (testing).         
%   tags = [cell of strings, Optional] tags of files saved for each condition. 
%          If not specified, infer from .coverage.txt.gz files in filedir
%
% Outputs
%   all_counts_3d = [cell of Nconditions arrays, size Nres x Nres x Nseq] for each tag,
%                   an Nres x Nres array counting the number of
%                   events (e.g., mutdel's) at each position given there is
%                   an event at the other position, for each of the
%                   reference sequences
%   all_counts    = [cell of Nconditions arrays, size Nseq x Nres] for each tag,
%                   the number of events (e.g., mutdel's) at each
%                   of the Nres positions, for each of the Nseq reference
%                   sequences. 
%   all_coverage  = [cell of Nconditions arrays, size Nseq x Nres] for each tag,
%                   the number of reads with coverage at each 
%                   of the Nres positions, for each of the Nseq reference
%                   sequences. 
%   tags = [cell of strings, Optional] tags of files read in.
%
% (C) Rhiju Das, HHMI & Stanford 2024-25

if ~exist( 'data_type', 'var')
    fprintf('Must specify data_type as "mutdel" (M2) or "del" (MOHCA)\n');
    return;
end

%% Read in all experimental conditions, just mutdel/mutdel
if ~exist('tags','var') | isempty(tags)
    files=dir([filedir,'/*.coverage.txt.gz']);
    %all_m2 = {}; all_counts = {}; all_coverage = {};
    tags = {};
    for i = 1:length(files)
        tags{i}=strrep(files(i).name,'.coverage.txt.gz','');
    end
end
for i = 1:length(tags)
    tag = tags{i};
    fprintf('Reading in info for tag: %s\n', tag )
    tic
    all_m2{i} = load_txt([filedir,'/',tag,'.counts_',data_type,'_',data_type,'.txt'],1);
    all_counts{i} = load_txt([filedir,'/',tag,'.counts_',data_type,'.txt'],1);
    all_coverage{i} = load_txt([filedir,'/',tag,'.coverage.txt'],1);
    toc
end

%% reshape
fprintf('Reshaping 2D data...\n');
Nseq = size(all_counts{1},1);
Nres = size(all_counts{1},2);
tic
for i = 1:length(tags)
    all_counts_3d{i} = reshape( all_m2{i}', Nres,Nres,Nseq);
end
toc
