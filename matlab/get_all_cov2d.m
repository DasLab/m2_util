function [all_cov2d, all_cov2d_err] = get_all_cov2d( all_counts_3d, all_counts, all_coverage, BLANK_OUT5, BLANK_OUT3 );
% [all_cov2d, all_cov2d_err] = get_all_cov2d( all_counts_3d, all_counts, all_coverage, BLANK_OUT5, BLANK_OUT3 );
%
% Derive cov2d = counts(i,j)/ counts(i) / counts(j) - 1
% 
% Inputs:
%
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
%  BLANK_OUT5 = Set to NaN this number of 5' residues. (Default 0)
%  BLANK_OUT3 = Set to NaN number of 3' residues. (Default 0)
%
% Outputs:
%   all_cov2d = [cell of Nconditions arrays, size Nres x Nres x Nseq] cov2d for
%                   each tag
%   all_cov2d_err = [cell of Nconditions arrays, size Nres x Nres x Nseq] error on cov2d for
%                   each tag
%
% (C) Rhiju Das, Stanford and HHMI, 2025

if ~exist('BLANK_OUT5','var') | isempty(BLANK_OUT5); BLANK_OUT5 = 0; end;
if ~exist('BLANK_OUT3','var') | isempty(BLANK_OUT3); BLANK_OUT3 = 0; end;

Nseq = size(all_counts{1},1);
Nres = size(all_counts{1},2);
which_res = [(1+BLANK_OUT5):(Nres-BLANK_OUT3)];
all_cov2d = {};
all_cov2d_err = {};
for n = 1:length(all_counts)
    r3=all_counts_3d{n};
    r1=all_counts{n};
    cov2d = NaN*ones(Nres,Nres,Nseq);
    cov2d_err = NaN*ones(Nres,Nres,Nseq);
    for idx = 1:Nseq
        cvg=max(all_coverage{n}(idx,:)); % would be better if coverage was a 2D matrix?
        c = cvg * r3(:,:,idx)./(r1(idx,:)'*r1(idx,:))-1;       
        cov2d(which_res,which_res,idx) = c(which_res,which_res);

        c_err = cvg * sqrt(1+r3(:,:,idx))./(r1(idx,:)'*r1(idx,:)); % assume noise in 2D counts dominate error
        cov2d_err(which_res,which_res,idx) = c_err(which_res,which_res);
    end
    all_cov2d{n} = cov2d;
    all_cov2d_err{n} = cov2d_err;
end