function [all_Zscores,all_Zscores_err] = get_all_Zscores( all_m2, all_m2_err );
% [all_Zscores,all_Zscores_err] = get_all_Zscores( all_m2, all_m2_err );
%
% Inputs:
%  all_m2     = [cell of Nconditions arrays, size Nres x Nres x Nseq] 
%                     chemical mapping profiles for each mutant
%  all_m2_err = [cell of Nconditions arrays, size Nres x Nres x Nseq] error
%                     on data to plot (Optional)
%
% Outputs
%
%  all_Zscores     = [cell of Nconditions arrays, size Nres x Nres x Nseq] 
%                     Zscores
%  all_Zscores_err = [cell of Nconditions arrays, size Nres x Nres x Nseq] 
%                     error on Zscores (set to nan if all_m2_err not
%                     provided)
%
% (C) R. Das, Stanford University, HHMI
if ~exist('all_m2_err','var') | isempty(all_m2_err)
    for n = 1:length(all_m2); all_m2_err{n} = nan*all_m2{n}; end;
end

all_Zscores = {};
all_Zscores_err = {};
for n = 1:length(all_m2)
    fprintf('Doing condition %d/%d\n',n,length(all_m2))
    m2 = all_m2{n};
    m2_err = all_m2_err{n};
    Zscores = 0*m2;
    Zscores_err = nan*m2;
    for j = 1:size(m2,3)
        [Zscores(:,:,j),Zscores_err(:,:,j)] = get_Zscore(m2(:,:,j),m2_err(:,:,j));
    end
    all_Zscores{n} = Zscores;
    all_Zscores_err{n} = Zscores_err;
end
