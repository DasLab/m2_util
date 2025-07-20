function signal_to_noise = get_map_signal_to_noise( all_maps, all_maps_err, MIN_SEQ_DISTANCE)
% signal_to_noise = get_map_signal_to_noise( all_maps, all_maps_err, MIN_SEQ_DISTANCE)
%
% Inputs
%  all_maps     = [cell of Nconditions arrays, size Nres x Nres x Nseq] data to plot
%  all_maps_err = [cell of Nconditions arrays, size Nres x Nres x Nseq] error on data to plot
%
% Optional inputs:
%  MIN_SEQ_DISTANCE = mask out residues pairs with sequence less than or
%                         equal to this number (Default 4)
%
% Output
%  signal_to_noise = [Nseq x Nconditions] signal_to_noise estimates
%
% (C) R. Das, Stanford/HHMI 2025

if ~exist('MIN_SEQ_DISTANCE','var') | isempty(MIN_SEQ_DISTANCE); MIN_SEQ_DISTANCE = 4; end;

for n = 1:length( all_maps )
    for j = 1:size( all_maps{n}, 3)
        map2d = all_maps{n}(:,:,j); 
        map2d_err = all_maps_err{n}(:,:,j); 
        [ii,jj]=meshgrid(1:size(map2d,1));
        map2d(abs(ii-jj) <= MIN_SEQ_DISTANCE) = nan; % mask residues close to each other
        signal_to_noise(j,n) = ubr_estimate_signal_to_noise_ratio( abs(map2d(:)), map2d_err(:) );
    end
end