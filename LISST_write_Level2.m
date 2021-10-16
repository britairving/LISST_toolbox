function LISST_write_Level2(cfg, data_proc, meta_proc, data_grid, meta_grid)
%FUNCTION LISST_write_Level2
%
%  Syntax:
%    [cfg, dat] = LISST_write_Level2(cfg,dat)
%  
%  Description: QC'ed data, format depends on cfg.write_format
%
%  Refereces:
%     https://seabass.gsfc.nasa.gov/wiki/Data_Submission
%     https://seabass.gsfc.nasa.gov/wiki/stdfields#Table%20of%20Field%20Names%20and%20Units
%
%  Notes:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
if nargin > 3
  include_gridded = 1;
else
  include_gridded = 0;
end
dbstop if error % debug mode

%% 1 | Create submit folder with specifics on inversion model and background files
cfg.path.dir_submit_doc = fullfile(cfg.path.dir_submit_L2,'documents');
if ~exist(cfg.path.dir_submit_doc,'dir')
  mkdir(cfg.path.dir_submit_doc) % store associated documents/zscat files
end
%% 2 | move matlab datafile, zscat files, and raw data files into appropriate folders
if ~exist(cfg.path.dir_submit_L0,'dir')
  mkdir(cfg.path.dir_submit_L0);
end
% only copy if directory is empty
l0_dir = dir(cfg.path.dir_submit_L0);
if isempty(l0_dir)
  fprintf('Copying raw data files to: %s\n',cfg.path.dir_submit_L0)
  unique_datfiles = unique(data_proc.datfile);
  for nr = 1:numel(unique_datfiles)
    copyfile(fullfile(cfg.path.dir_raw,unique_datfiles{nr}),cfg.path.dir_submit_L0)
  end
end
fprintf('Copying background files to: %s\n',cfg.path.dir_submit_doc)
unique_zscfiles = unique(data_proc.zscfile);
for nz = 1:numel(unique_zscfiles)
  copyfile(fullfile(cfg.path.dir_zsc,unique_zscfiles{nz}),cfg.path.dir_submit_doc)
end

% fprintf('Copying MATLAB files to: %s\n',cfg.path.dir_submit_doc)
% copyfile(cfg.path.file_qc,cfg.path.dir_submit_doc)
% if include_gridded
%   copyfile(cfg.path.file_grid,cfg.path.dir_submit_doc)
% end
% 

%% 2 | Write processed data to 
switch cfg.write_format
  %% I  | SeaBASS format
  case 'seabass' % https://seabass.gsfc.nasa.gov/wiki/Data_Submission
    LISST_write_SeaBASS(cfg,data_proc,meta_proc,data_grid)
  %% II | NGA LTER Format
  case 'nga_lter' % AXIOM Research Workspace
    LISST_write_NGALTER(cfg,data_proc,meta_proc)
  case 'netcdf'
    fprintf('This write format is set up yet\n')
    keyboard
  otherwise
    fprintf('This write format is set up yet\n')
    keyboard
end

end %% MAIN FUNCTION 