function  m = load_txt( filename, verbose);
% m = load_txt( filename );
% 
% Read text file with counts as integers, can be comma-separated, and can
%   be gzipped.
%
% Note: This is a clone of READ_TABLE_FILE from UBR.
%
% Input
%   filename = txt or txt.gz with counts as integers
%   verbose  = print out filename (Default 0)
%
% Output
%   m = table in matlab format.
%
% (C) R. Das, Stanford University, HHMI 2023-25
 
if ~exist('verbose','var') verbose = 0; end;

gzip_file = [filename,'.gz'];
if ~exist(filename,'file')
    %fprintf( 'Unzipping %s... \n',gzip_file);
    if ~exist(gzip_file,'file') fprintf( '\n\nDid not find %s or %s! \n\n',filename,gzip_file); end;
    assert(exist(gzip_file,'file'));
    gunzip( gzip_file );
end

if verbose; fprintf('Reading in %s\n',filename); end

t = readtable(filename);
m = table2array(t,'unt32');

if exist(gzip_file,'file')
    delete(filename); % to save space.
end
