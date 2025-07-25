function coverage = get_coverage( all_coverage, SHOW_PLOT )
% coverage = get_coverage( all_coverage )
%
% Inputs:
% all_coverage  = [cell of Nconditions arrays, size Nseq x Nres] for each tag,
%                   the number of reads with coverage at each 
%                   of the Nres positions, for each of the Nseq reference
%                   sequences, as saved by bam_to_m2.
% SHOW_PLOT = make a plot of coverage vs. position for first condition.
%
% Outputs:
% coverage = [Nseqs x Nconditions array] coverage for each sequence across each
%                 of the conditions
%
%
if ~exist('SHOW_PLOT','var'); SHOW_PLOT = 0; end;
coverage = [];
for n = 1:length(all_coverage); 
    coverage(:,n) = max(all_coverage{n},[],2);
end

if SHOW_PLOT
    figure(3);
    set(gcf,'color','white'); clf
    [~,sortidx] = sort(sum(coverage'),'descend');
    imagesc(  all_coverage{1}(sortidx',:),[0 max( all_coverage{1}(:))])
    xlabel('Position');
    ylabel('Sequence (ordered by coverage)');
    title( 'Coverage')
    colormap(parula)
    colorbar();
end
