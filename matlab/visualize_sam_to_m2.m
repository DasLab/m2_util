function  visualize_sam_to_m2( tag, maxscale, seqpos, SAVE_FIGURE, sequence_names );
%  visualize_sam_to_m2( tag, maxscale, seqpos );
%  visualize_sam_to_m2( dir, maxscale, seqpos );
%
% Inputs
%  tag = (string) output tag of files from sam_to_m2.py. Files like
%        tag.counts_del_del.txt, ... should exist
%     OR
%  dir = (string) directory name where to look for files like
%         *counts_del_del.txt.  [Default: current directory]
%  maxscale = maximum value to show in map [default 0.05]
%  seqpos = conventional numbering to apply to sequences.
%  SAVE_FIGURE = save figures to Figures/ directory
%
% (C) R. Das, HHMI & Stanford University 2025

if ~exist('tag','var') | isempty(tag), tag = './'; end; % look in current directory
if ~exist('maxscale','var') | maxscale==0; maxscale = 0.05; end;
if ~exist('seqpos','var'); seqpos = []; end;
if ~exist('SAVE_FIGURE','var'); SAVE_FIGURE = 1; end;
if ~exist('sequence_names','var'); sequence_names = {}; end;

if exist(tag,'dir') & ~exist( [tag,'.counts_del.txt'], 'file')
    data_dir=tag;
    files = dir([data_dir,'*.counts_del.txt']);
    for n = 1:length(files); cols = strsplit(files(n).name,'.'); tags{n} = cols{1}; end
    for n = 1:length(tags)
        visualize_sam_to_m2( tags{n}, maxscale, seqpos, SAVE_FIGURE, sequences );
    end
    return;
end

%% Mutate-and-map?
m2_del_del_all = load_txt([tag,'.counts_del_del.txt']);
m2_del_mut_all = load_txt([tag,'.counts_del_mut.txt']);
m2_mut_del_all = load_txt([tag,'.counts_mut_del.txt']);
m2_mut_mut_all = load_txt([tag,'.counts_mut_mut.txt']);
del_all = load_txt([tag,'.counts_del.txt'])';
mut_all = load_txt([tag,'.counts_mut.txt'])';
cov_all = load_txt([tag,'.coverage.txt']);
Nres = size(m2_del_del_all,2);
Nseq = size(m2_del_del_all,1)/Nres;
fprintf('\n');
for i = 1:Nseq;
    idx = Nres*(i-1)+ [1:Nres];
    m2_del_del = m2_del_del_all(idx,:);
    m2_del_mut = m2_del_mut_all(idx,:);
    m2_mut_del = m2_mut_del_all(idx,:);
    m2_mut_mut = m2_mut_mut_all(idx,:);
    cov = cov_all(i,:);
    del = del_all(:,i);
    mut = mut_all(:,i);
    sequence_name = '';
    if i <= length(sequence_names);  sequence_name=sequence_names{i}; end;
    fprintf('%8d deletions (%7.4f) and %8d mutations for (%7.4f) coverage %8d  %s %s\n', sum( del ), sum(del)/max(cov), sum(mut), sum(mut)/max(cov), max(cov), tag, sequence_name);
    if isempty(seqpos); seqpos = [1:length(del)]; end;

    set(figure(4),'color','white','pos', [50 50 1000 900]);
    clf;
    subplot(2,2,1);
    imagesc(seqpos,seqpos,m2_del_del./del,[0 maxscale]);
    colormap(1-gray(100)); xlabel('Mod position'); ylabel('Anchor position'); colorbar();
    title({'Delete and delete',sequence_name,tag},'Interpreter','none')

    subplot(2,2,2);
    imagesc(seqpos,seqpos,m2_del_mut./del,[0 maxscale]);
    colormap(1-gray(100)); xlabel('Mod position'); ylabel('Anchor position'); colorbar();
    title({'Delete and Mutate',sequence_name,tag},'Interpreter','none')

    subplot(2,2,3);
    imagesc(seqpos,seqpos,m2_mut_del./mut,[0 maxscale]);
    colormap(1-gray(100)); xlabel('Mod position'); ylabel('Anchor position'); colorbar();
    title({'Mutate and Delete',sequence_name,tag},'Interpreter','none')

    subplot(2,2,4);
    imagesc(seqpos,seqpos,m2_mut_mut./mut,[0 maxscale]);
    colormap(1-gray(100)); xlabel('Mod position'); ylabel('Anchor position'); colorbar();
    title({'Mutate and Mutate',sequence_name,tag},'Interpreter','none')

    if SAVE_FIGURE
        outdir = 'Figures/';
        if ~exist(outdir,'dir'); mkdir( outdir ); end;
        [~,out_tag] = fileparts(tag);
        outfile = ['Figures/',out_tag,'.png'];
        if length(sequence_name)>0
            outfile = ['Figures/',out_tag,'_',sequence_name,'.png'];
        end
        fprintf('Creating: %s\n',outfile);
        figure(4); print( outfile, '-dpng');
    end

    %%
    % set(figure(5),'color','white');
    % clf
    % m2_mutdel_mutdel = load([tag,'.counts_mutdel_mutdel.txt']);
    % mutdel = load([tag,'.counts_mutdel.txt'])';
    % imagesc(seqpos,seqpos,m2_mutdel_mutdel./mutdel,[0 maxscale]);
    % colormap(1-gray(100)); xlabel('Map position'); ylabel('Mutate position')
    % title({'Mut/del and mut/del',tag},'Interpreter','none')
    %
    
    %if ~SAVE_FIGURE & (i < Nseq); pause; end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = load_txt(infile);
% load infile -- gunzip if necessary.
x = [];
infile_gz = '';
if ~exist(infile,'file')
    infile_gz = [infile,'.gz'];
    assert( exist( infile_gz,'file' ));
    %system(['gunzip -k ',infile_gz] );
    gunzip(infile_gz)
    assert( exist( infile,'file' ));
end
x = load( infile );
if length(infile_gz)>0
    delete(infile);
end
