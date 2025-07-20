function output_maps_to_hdf5( all_maps, dataset_name, BLANK_OUT5, BLANK_OUT3, reads, sequences, headers, conditions, idx, which_conditions, all_maps_err, signal_to_noise)
% output_maps_to_hdf5( all_maps, dataset_name, BLANK_OUT5, BLANK_OUT3,reads, sequences, headers, conditions, idx, which_conditions)
%
% Inputs
% all_maps     = [cell of Nconditions arrays, size Nres x Nres x Nseq] data to plot
% dataset_name = name of outfile (suffix will be added: .2D.hdf5')
% BLANK_OUT5   = Number of 5' residues missing from data
% BLANK_OUT3   = Number of 3' residues missing from data 
% reads     = [Nseq x Nconditions] total reads contributing to each map (for
% sequences    = [cell of Nseq strings] sequences of A,C,G,U
% headers      = [cell of Nseqs strings] short name of sequences for output to map
% conditions   = [cell of Nconditions strings] short name of each experimental condition
%
% Optional inputs
% idx             = [list of integers] indices of which sequences to plot. Default:
%                  all sequences
% which_conditions= [list of integers] indices of which conditions to show.
% all_maps_err    = [Nseq x Nconditions] signal_to_noise estimates. 
%                      Default []
% signal_to_noise = [Nseq x Nconditions] signal_to_noise estimates. 
%                      Default []
%
% (C) R. Das, Stanford/HHMI 2025

if ~exist('idx','var') | isempty(idx); idx = [1:length(sequences)]; end;
if ~exist('which_conditions','var') | isempty(which_conditions); which_conditions = [1:length(all_maps)]; end;
if ~exist('all_maps_err','var') | isempty(all_maps_err)
    for n = 1:length(all_maps); all_maps_err{n} = nan*all_maps{n}; end;
end
if ~exist('signal_to_noise','var'); signal_to_noise = []; end;

assert( length(conditions) == length( all_maps ));

d = struct();
d.dataset_name=dataset_name;
d.sequences = sequences(idx);
d.headers = headers(idx);
d.BLANK_OUT5 = BLANK_OUT5;
d.BLANK_OUT3 = BLANK_OUT3;
d.conditions = conditions(which_conditions);
d.reads = reads(idx,which_conditions);
%d.total_reads = sum(reads(:,which_conditions),2);
d.signal_to_noise = signal_to_noise(idx,which_conditions); 
tic
for q = 1:length(which_conditions)
    fprintf('Reshaping %s\n',conditions{which_conditions(q)});
    d.r_norm(:,:,:,q) = permute( all_maps{which_conditions(q)}(:,:,idx), [3,1,2] );
    d.r_norm_err(:,:,:,q) = permute( all_maps_err{which_conditions(q)}(:,:,idx), [3,1,2] );
end
toc
output_ubr_hdf5( [dataset_name,'.2D.hdf5'], d, 1);
