function  [Zscore, Zscore_err] = get_Zscore(m2, m2_err);
% [Zscore, Zscore_err] = get_Zscore(m2, m2_err);
%
% m2 = [Nres x Nres] chemical mapping profile for mutation at each position
% m2_err = (Optional) [Nres x Nres] error estimate
%
% (C) R. Das, Stanford University & HHMI, 2025

if ~exist('m2_err','var') | isempty(m2_err); m2_err = nan*m2; end;

Zscore = [];
Zscore_err = [];
for i = 1:size(m2,1)
    s = nanstd( m2(i,:) );
    m = nanmean( m2(i,:) );
    Zscore(i,:) = (m2(i,:) - m)/s;
    Zscore_err(i,:) = m2_err(i,:)/s;
end
