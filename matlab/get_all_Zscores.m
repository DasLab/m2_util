function all_Zscores = get_all_Zscores( all_m2 );
% all_Zscores = get_all_Zscores( all_m2 );

all_Zscores = {};
for n = 1:length(all_m2)
    fprintf('Doing condition %d/%d\n',n,length(all_m2))
    m2 = all_m2{n};
    Zscores = 0*m2;
    for j = 1:size(m2,3)
        Zscores(:,:,j) = get_Zscore(m2(:,:,j));
    end
    all_Zscores{n} = Zscores;
end
