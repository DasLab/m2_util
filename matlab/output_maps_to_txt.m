function output_maps_to_txt(all_maps,tags,outdir);
% output_maps_to_txt(all_maps,tags,outdir);
%
% Output M2 maps (e.g., cov2d) out to .2D.txt files.
%
% Note: OUTPUT_MAPS_TO_HDF5 is preferred...
%
% Inputs
%   all_maps = [cell of Nconditions arrays, size Nres x Nres x Nseq] data to plot
%       tags = [cell of Nconditions strings] short name of each experimental condition
%
% Optional input:
%   outdir  = output director (Default is './')
%

if ~exist('outdir','var') | isempty(outdir); outdir = './'; end;

Nres = size(all_maps{1},1);
if ~exist(outdir,'dir'); mkdir(outdir); end;
for n = 1:length(tags)
    out = reshape(all_maps{n},Nres,[])';
    outfile = [outdir,'/',tags{n},'.2D.txt'];
    fprintf('Outputting file: %s\n',outfile);
    t=table(out);
    writetable(t,outfile,'WriteVariableNames',false);
end
