function LISST_write_Level1(cfg, data_raw, meta_raw)
%FUNCTION LISST_write_Level1
%
%  Syntax:
%    [cfg, dat] = LISST_write_Level1(cfg,dat)
%  
%  Description: Raw instrument data, ASCII format
%
%  Notes:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Generate filename for writing
switch cfg.write_format
  % I | NGA LTER Format
  case 'nga_lter' % AXIOM Research Workspace
    wname = ['NGA_' cfg.header.cruise_ID '_LISST_L1_v1.csv'];
    wname = fullfile(cfg.path.dir_submit_L1,wname);
  otherwise
    wname = fullfile(cfg.path.dir_submit_L1,[cfg.project '_L1.csv']);
end

%% Generate appropriate fieldnames from column descriptions
cols = strrep(meta_raw.data.fields,'*','x');
cols = strrep(cols,'+','_');

%% Convert to table
raw = array2table(data_raw.data);
raw.Properties.VariableNames = cols;
raw.datfile = data_raw.datfile;

%% Format header for writing
hdr = strjoin(raw.Properties.VariableNames,',');

%% Specify file writing format
fmt = repmat( '%d,', 1,40);
fmt = [fmt '%s\n'];

% reformat data table to temporary variable to enable simply write to file
raw_write = table2cell(raw); % convert to cell array to handle different variable types
raw_write = raw_write';       % transpose because fprintf function prints data columnwise
% % convert NaNs to missing identifier
% raw_write(cellfun(@(x) isnumeric(x) && isnan(x), raw_write)) = {-9999}; % missing=-9999 SeaBASS required header field

%% Write odv table to NONDIFFERENTIAL file
fprintf('\n  Writing LISST data to: %s\n',wname);
fileID = fopen(wname,'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',wname)
  keyboard
end
fprintf(fileID,'%s\n',hdr); % write header
fprintf(fileID,fmt,raw_write{:});      % write data
fclose(fileID);                 % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
%fprintf('%s\n',cols{:}) % write header
%fprintf(fmt,raw_write{:})      % write data



end %% MAIN FUNCTION 