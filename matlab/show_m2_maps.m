function show_m2_maps( all_maps, idx, colorscale, BLANK_OUT5, BLANK_OUT3, coverage, headers, tags, descriptions, which_tags )
% show_m2_maps( all_maps, idx, colorscale, BLANK_OUT5, BLANK_OUT3, coverage, headers, tags, descriptions, which_tags )
%
% Show M2 maps, with different experimental conditions tiled, and each 
%   sequence in a new plot.  Pauses between different sequences.
%
% Inputs
%   all_maps = [cell of Nconditions arrays, size Nres x Nres x Nseq] data to plot
%
% Optional inputs
%  idx    = [list of integers] indices of which sequences to plot
%  colorscale  = maxscale, or [-minscale, maxscale], or cell of such
%                   numbers for each of Nconditions condition. Default: [0,0.05].
%  BLANK_OUT5 = Set axes to crop out this number of 5' residues. (Default 0)
%  BLANK_OUT3 = Set axes to crop out this number of 3' residues. (Default 0)
%  coverage = [Nseq x Nconditions] total reads contributing to each map (for
%     output in map)
%  headers = [cell of Nseqs strings] short name of sequences for output to map
%  tags = [cell of Nconditions strings] short name of each experimental condition
%  descriptions = [cell of Nseqs strings] longer description of sequences for output to
%              command window
%  which_tags = [list of integers] indices of which conditions to show.
%
% (C) R. Das, Stanford & HHMI, 2025

Nconditions = length( all_maps );
if Nconditions == 0; return; end;
Nres = size(all_maps{1},1);
Nseq = size(all_maps{1},3);

% how to tile plots; adjust as needed.
if Nconditions == 1; n1 = 1; n2 = 1; pix1 = 500; pix2 = 500;
elseif Nconditions == 2; n1 = 1; n2 = 2;  pix1 = 1000; pix2 = 500;
elseif Nconditions == 3; n1 = 1; n2 = 3;  pix1 = 1500; pix2 = 500;
elseif Nconditions == 4; n1 = 2; n2 = 2;  pix1 = 1000; pix2 = 1000;
else n1 = ceil(Nconditions/2); n2 = 2; pix1 = 1000; pix2 = 5000;
end

if ~exist('idx', 'var') | isempty(idx); idx = [1:Nseq]; end;
if size(idx,2)==1; idx = idx'; end;

if ~exist('colorscale', 'var' ) | isempty(colorscale) | colorscale==0; 
    colorscale = [0 0.05]; 
    warning('colorscale not specified, so using [0,0.05]. If plotting cov2d data, suggest using [-5,5].\n')
end;
if isnumeric(colorscale) & length(colorscale)==1; colorscale = [0, colorscale]; end;
if ~iscell( colorscale), colorscale  = {colorscale}; end;
if length(colorscale) == 1; colorscale = repmat(colorscale,1,Nconditions); end;

if ~exist('BLANK_OUT5','var') | isempty(BLANK_OUT5); BLANK_OUT5 = 0; end;
if ~exist('BLANK_OUT3','var') | isempty(BLANK_OUT3); BLANK_OUT3 = 0; end;
if ~exist( 'coverage', 'var') | isempty(coverage); coverage = nan * ones(Nseq,Nconditions); end;
if ~exist( 'headers', 'var')
    for n = 1:Nseq; headers{n} = sprintf('Sequence %d',n); end;
end;
if ~exist( 'descriptions', 'var')
    for n = 1:Nseq; descriptions{n} = ''; end;
end;
if ~exist( 'tags', 'var') | isempty( tags )
    for n = 1:Nconditions; tags{n} = sprintf('Condition %d',n); end;
end
if ~exist( 'which_tags', 'var') | isempty(which_tags), which_tags = [1:length(tags)]; end;

set(gcf,'color','white','Position',[100 100 pix1 pix2])
axislim = [1+BLANK_OUT5, Nres-BLANK_OUT3];
for i = 1:length(idx)
    k = idx(i);
    for q = which_tags
        n = which_tags(q);
        subplot(n1,n2,n);

        imagesc(all_maps{n}(:,:,k)',colorscale{n});

        plot_title = {...
            ['Sequence ',num2str(k)],...
            headers{k}(1:min(50,end)),...
            tags{n}};
        if ~isnan(coverage(k,n)); plot_title{1} = ['Sequence ',num2str(k), ' Counts ',num2str(coverage(k,n))]; end;

        title(plot_title,'Interpreter','none','fontsize',14);
        xlim( axislim ); ylim( axislim );
    end

    if colorscale{n}(1) < 0; colormap( [0.8 0.8 0.8; redwhiteblue(-1,1) ]);
    else colormap(1-gray(100)); end;

    fprintf('\nEntry %d/%d, index: %d, total counts: %d\n%s\n%s\n',i,length(idx),k,sum(coverage(k,:)),headers{k},descriptions{k})

    if k ~= idx(end); pause; end;
end

