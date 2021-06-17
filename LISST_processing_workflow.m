function LISST_processing_workflow
%% function LISST_processing_workflow
%  Description:
%    Wrapper script for processing LISST raw (binary .DAT) files.
%    Detailed documentation found in each script.
%
%  REQUIRED INPUT:  Need to change cfg.project and cfg.year.
%                   Organize data following LISST_process_options documentation
%                   
%  Refereces:
%    http://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/
%    LISST-Deep-Users-Manual-May-2013.pdf
%    http://www.sequoiasci.com/article/how-lisst-instruments-measure-the-size-distribution-and-concentration-of-particles/
%
%  Notes:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%   | Starting point


%% 0 | USER INPUT: project and year
% cfg.project      = 'LISST_sn4025_2018_EXPORTS_SR201812';
cfg.project      = 'LISST_sn4041_2018_ASGARD_SKQ201813S';     % Project folder in LISST_Data
% cfg.project      = 'LISST_sn4025_2017_ASGARD_SKQ201709S';     % Project folder in LISST_Data
cfg.year         = 2018;                                      % year when first measurement was taken
cfg.testing      = 0;                                         % 0 = processes all available DAT files, 1 = processes first 10 DAT files
cfg.savefig      = 1;                                         % 0 = does not save figures, 1 = saves figures to cfg.path.dir_figs
cfg.path.base    = '/Volumes/toshiba1/LISST_Data_Structured'; % Path to where LISST_Data folder
cfg.path.toolbox = '/Users/bkirving/Documents/MATLAB/LISST_toolbox/'; % Path to LISST toolbox

%% 0 | Configure processing
skip_to_proc =  0;    % 1 = jumps to processing
skip_to_qaqc =  0;    % 1 = jumps to qaqc
%% 1 | Read project metadata and instrument information
try % Try to read in project configuration from [cfg.project '_config.m]
  addpath(genpath(cfg.path.base))
  eval(['cfg = ' cfg.project '_config(cfg);'])
catch % catch and explain why script stopped
  fprintf('Need to set up %s_config.m script, see example provided\n',cfg.project)
  error('Cannot load config.m for this project\n')
end

%% 2 | Configure paths and processing methods
cfg = LISST_processing_config(cfg);
%% 3 | step through processing workflow
if ~skip_to_qaqc
  if ~skip_to_proc
    %% 4 | Read raw data
    % Reads raw data into structure
    [cfg, data_raw, meta_raw] = LISST_read_raw_data(cfg);

    %% 5 | Write raw data to ASCII file
   % if strcmpi(cfg.write_format,'nga_lter') 
   %   LISST_write_Level1(cfg,data_raw,meta_raw); % Raw instrument data, ASCII format
   % end
   
    %% 6 | Preprocess data
    % Identifies downcasts & matched with CTD casts
    [cfg, data_pre, meta_pre] = LISST_preprocess_data(cfg, data_raw, meta_raw);
  else
    fprintf('Loading preprocessed data from file: %s\n',cfg.path.file_pre)
    load(cfg.path.file_pre);            % load raw data in Matlab format
    cfg.path.toolbox = '/Users/bkirving/Documents/MATLAB/LISST_toolbox/'; % Path to LISST toolbox
    cfg = LISST_processing_config(cfg); % reload load cfg incase things have changed
  end
  
  %% 7 | Process data
  % Processes data following Sequoia's recommendations and derives various other parameters. 
  [cfg, data_proc, meta_proc] = LISST_process_data(cfg,data_pre,meta_pre);
  
else
  fprintf('Loading processed data from file: %s\n',cfg.path.file_proc)
  load(cfg.path.file_proc);           % load preprocessed data 
  cfg = LISST_processing_config(cfg); % reload load cfg incase things have changed
end

%% 8 | Automated QC LISST data
% Flags data according to QC tests described in LISST manual. 
[cfg, data_proc, meta_proc] = LISST_data_qaqc_auto(cfg, data_proc, meta_proc);
  
%% 9 | Grid data
% Grids data 
[cfg, data_grid, meta_grid] = LISST_grid_data(cfg,data_proc,meta_proc);

%% 10 | Manual QC LISST data
% NOT CREATED YET! Will.. Loop through casts and flag data manually. 
if cfg.qaqc_options.expert_qc
  fprintf('SET UP EXPERT QC!\n')
  keyboard
  [cfg, data_proc, meta_proc] = LISST_data_qaqc_expert(cfg, data_proc, meta_proc);
end

%% 11 | Write to file
% Organizes depend files into folders, and writes processed data to text
% file depending on cfg.write_fmt
LISST_write_Level2(cfg, data_proc, meta_proc, data_grid, meta_grid);
keyboard

%% 12 | Visualize
% Plots gridded data
% LISST_plot_data_grid(cfg,dat,dat_bin,data);
