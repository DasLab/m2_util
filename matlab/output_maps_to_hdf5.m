function output_maps_to_hdf5( all_maps, dataset_name, BLANK_OUT5, BLANK_OUT3, coverage, sequences, headers, tags, idx, which_tags)
% output_maps_to_hdf5( all_maps, dataset_name, BLANK_OUT5, BLANK_OUT3,coverage, sequences, headers, tags, idx, which_tags)
%
% Inputs
% all_maps     = [cell of Nconditions arrays, size Nres x Nres x Nseq] data to plot
% dataset_name = name of outfile (suffix will be added: .2D.hdf5')
% BLANK_OUT5 = Number of 5' residues missing from data
% BLANK_OUT3 = Number of 3' residues missing from data 
% coverage = [Nseq x Nconditions] total reads contributing to each map (for
% sequences    = [cell of Nseq strings] sequences of A,C,G,U
% headers = [cell of Nseqs strings] short name of sequences for output to map
% tags         = [cell of Nconditions strings] short name of each experimental condition
%
% Optional inputs
% idx    = [list of integers] indices of which sequences to plot. Default:
%               all sequences
% which_tags = [list of integers] indices of which conditions to show.
%
% (C) R. Das, Stanford/HHMI 2025

if ~exist('idx','var') | isempty(idx); idx = [1:length(sequences)]; end;
if ~exist('which_tags','var') | isempty(which_tags); which_tags = [1:length(all_maps)]; end;

assert( length(tags) == length( all_maps ));

d = struct();
d.dataset_name=dataset_name;
d.sequences = sequences(idx);
d.headers = headers(idx);
d.BLANK_OUT5 = BLANK_OUT5;
d.BLANK_OUT3 = BLANK_OUT3;
d.tags = tags(which_tags);
d.coverage = coverage(idx,which_tags);
d.total_coverage = sum(coverage(:,which_tags),2);
tic
for q = 1:length(which_tags)
    fprintf('Reshaping %s\n',tags{which_tags(q)});
    d.cov2d(:,:,:,q) = permute( all_maps{which_tags(q)}(:,:,idx), [3,1,2] );
end
toc
output_ubr_hdf5( [dataset_name,'.2D.hdf5'], d, 1);
