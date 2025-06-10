function  Zscore = get_Zscore(m2);

Zscore = [];
for i = 1:size(m2,1)
    s = nanstd( m2(i,:) );
    m = nanmean( m2(i,:) );
    Zscore(i,:) = (m2(i,:) - m)/s;
end
