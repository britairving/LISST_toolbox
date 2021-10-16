function cfg = LISST_processing_config(cfg)
%function LISST_processing_config
%
%  Syntax:
%    cfg = LISST_processing_config(cfg)
%
%  Description:
%   ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
%          Go through this script carefully to set up processing 
%   ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **
%
%    Builds paths to LISST data and instrument files based on expected
%    folder hierarchy.
%    Reads/defines instrument information important for data processing and
%    metadata storage. 
%    Aside from path_length, and inversion model used to process data,
%    variables are read from instrument files lisst.ini and
%    InstrumentData.txt using paths built in the beginning part of this
%    script.
%    Defines QC flags and defines automatic QC tests. 
%
%    Expected folder hierarchy...
%     project  = Path to project folder
%     dir_raw  = Path to project binary data files (*.DAT)
%     dir_proc = Path to folder where all matlab files are saved
%     dir_zsc  = Path to project background scattering files
%     dir_inst = Path to project instrument information (Factory_zsc_[S/N].asc, InstrumentData.txt, lisst.ini, Ringarea_[S/N].asc)
%     dir_figs = Path to folder where figures will save to
%     dir_ctd  = Path to CTD data files, formatting defined by cfg.ctd_type
%     ---------------------------------------------------------------------
%     Example of expected folder hierarchy: organize data as follows
%       cfg.path.base\
%         LISST_sn4025_2018_exports_np_sr1812\ | (project)
%           background_scatter_files\          | (dir_zsc)
%           instrument_info\                   | (dir_inst)
%           raw\                               | (dir_raw)
%           proc\                              | (dir_proc)
%           figures\                           | (dir_figs)
%           ctd_files\                         | (dir_ctd)
%           submit\L0\                         | (dir_submit_L0)
%           submit\L1\                         | (dir_submit_L1)
%           submit\L2\                         | (dir_submit_L2)
%     ---------------------------------------------------------------------

%  References:
%    http://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/
%    LISST-Deep-Users-Manual-May-2013.pdf
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
% -----------------------------------------------------------------------
%% 1 | Configure paths and files 
% -----------------------------------------------------------------------
% This is all build on cfg.path.base, defined in
% LISST_processing_workflow.m
cfg.path.project  = fullfile(cfg.path.base,cfg.project);                   % Path to project folder
cfg.path.dir_raw  = fullfile(cfg.path.project,'raw');                      % Path to project binary data files (*.DAT)
cfg.path.dir_zsc  = fullfile(cfg.path.project,'background_scatter_files'); % Path to project background scattering files
cfg.path.dir_inst = fullfile(cfg.path.project,'instrument_info');          % Path to project instrument information (Factory_zsc_[S/N].asc, InstrumentData.txt, lisst.ini, Ringarea_[S/N].asc, LISST-100X(and_LISST-100)_spreadsheet_for_date_conversion.xls)
cfg.path.dir_figs = fullfile(cfg.path.project,'figures');                  % Path to project figures
cfg.path.dir_ctd  = fullfile(cfg.path.project,'ctd_files');                % Path to project ctd data files (UNBINNED or Level0/1! Need downcast and upcast for easy cast matching)
cfg.path.dir_proc = fullfile(cfg.path.project,'proc');                     % Path to processed data files
cfg.path.dir_submit_L0 = fullfile(cfg.path.project,'submit','L0');         % Path to raw data in original format
cfg.path.dir_submit_L1 = fullfile(cfg.path.project,'submit','L1');         % Path to raw data in ascii format
cfg.path.dir_submit_L2 = fullfile(cfg.path.project,'submit','L2');         % Path to processed qc'd data in ascii format

%%  1a | Update paths and files based on the above specifications
cfg.path.file_ringarea = fullfile(cfg.path.dir_inst,['Ringarea_' num2str(cfg.inst.LISSTsn)    '.asc']);  % Path to ring area file
cfg.path.file_factoryz = fullfile(cfg.path.dir_inst,['Factory_zsc_' num2str(cfg.inst.LISSTsn) '.asc']);  % Path to factory zscat (background scatter) file
cfg.path.file_ctddata  = fullfile(cfg.path.dir_ctd, ['CTD_' cfg.ctd_type '.mat']);                       % Path to file that contains CTD data in matlab format
cfg.path.file_ctdmatch = fullfile(cfg.path.dir_proc,[cfg.project '_CTD_cast_match.csv']);                % Path to file that contains CTD cast number matched with LISST profile
cfg.path.file_ctdlag   = fullfile(cfg.path.dir_proc,[cfg.project '_CTD_lag_corrections.csv']);           % Path to file that contains time and depth lags compared to CTD data
cfg.path.file_downcast = fullfile(cfg.path.dir_proc,[cfg.project '_downcast_indices.csv']);              % Path to file that contains indicies for start and end of downcast for cast

cfg.path.file_raw      = fullfile(cfg.path.dir_proc,[cfg.project '_raw.mat']);  % STEP 1: Raw data read into matlab strucutre
cfg.path.file_pre      = fullfile(cfg.path.dir_proc,[cfg.project '_pre.mat']);  % STEP 2: Preprocessed data - matched with CTD casts, downcasts identified, temperature and depth corrected for any lags

%%  1b | Create folder hierarchy if doesn't exist
fields = {'dir_raw' 'dir_zsc' 'dir_inst' 'dir_figs' 'dir_ctd' 'dir_proc'};
for nf = 1:numel(fields)
  if ~exist(cfg.path.(fields{nf}),'dir')
    fprintf(' making directory %s\n',cfg.path.(fields{nf}))
    fprintf(' ** need to put appropriate files into directory %s!\n',cfg.path.(fields{nf}))
    keyboard
    mkdir(cfg.path.(fields{nf}))
  end
end

%%  1c | Add paths
addpath(genpath(cfg.path.toolbox));
addpath(fullfile(cfg.path.toolbox,'sequoia'))           % path to sequoia scripts (downloaded Feb 2020 from http://www.sequoiasci.com/wp-content/uploads/2015/01/matlab-processing.zip)
addpath(fullfile(cfg.path.toolbox,'lisst_utility'))     % path to random dependent files
if ~exist(cfg.path.dir_figs,'dir')
  mkdir(cfg.path.dir_figs)
end

% -----------------------------------------------------------------------
%% 2 | Processing methods based on user preference
% -----------------------------------------------------------------------
cfg.proc_options.header = 'structure to store all switches and settings used for processing LISST data';
cfg.proc_options.calculate_range = 0;   % 0 = does nothing, 1 = calculates range of parameters by taking the difference of the parameter derived using the maximum and minimum backgrounds
cfg.proc_options.profile_limit   = 7;   % [meter] ideally just below soak depth - depth profile limit for cast identification.  I.e. if data collected on shallow shelf, set low. 

%%  2b | Define inversion model used to process data
% The Spherical particle model performs the mathematical inversion under
% the assumption that the particles that scattered light are all spheres.
% Sequoia employs the so-called full Mie solution when using this approach.
% The Mie solution is a generalized solution to the scattering of light
% from spheres and is being used as the standard inversion method by all
% laser manufacturers.
% However, sediment particles in the aquatic environment are never perfect
% spheres. Consequently Sequoia provides an option to invert the scattering
% pattern under the assumption that the particles are randomly (or
% irregularly) shaped particles. The exact details of how this inversion
% takes place is described in a scientific paper by Agrawal et al.6, which
% can be downloaded from the library section on Sequoia’s website
% (www.SequoiaSci.com). The direct URL is
% https://www.sequoiasci.com/library/technical-papers/. 
% A brief version of the method is described in this News article on
% Sequoia’s website: http://sequoiasci.com/Articles/ArticlePage.aspx?pageI
% d=133
% So what method should you use? If you are going to compare with the
% output from another laser particle sizer, for example a laboratory
% particle sizer, you should choose the Spherical particle model. The
% reason for this is that no other laser manufacturer provides a randomly
% shaped particle model for inversion, so there will be nothing to compare
% to. Generally, it is recommended that you keep both boxes checked (this
% is also the default option), so that you won’t have to reprocess your
% data because you need to see what the data look like when processed as
% randomly shaped particles.
% RANDOM = 1 if matrices based on scattering from randomly shaped particles are to be used for inversion. NOTE: Only type B and C instruments are supported for RANDOM = 1.
% RANDOM = 0 if matrices based on scattering from spherical particles are to be used for inversion.
cfg.proc_options.rand = 0;

% Define processed filename based on model
if cfg.proc_options.rand == 1 % non-spherical inversion model
  model_str = 'nonspherical';
elseif cfg.proc_options.rand == 0 % spherical inversion model
  model_str = 'spherical';
end
cfg.path.file_proc     = fullfile(cfg.path.dir_proc,[cfg.project '_proc_' model_str '.mat']); % STEP 3: processed data - defined in 

% -----------------------------------------------------------------------
%% 3 | Gridding methods
% -----------------------------------------------------------------------
cfg.grid_options.bin_depth_m     = 1.0; % [meter] depth range for binning
cfg.grid_options.ignore_flagged  = 0;   % 1 = throws out flagged data before gridding, 0 = grids all data regardless if flagged or not

% -----------------------------------------------------------------------
%% 4 | Processing methods based on LISST instrument specifications
% -----------------------------------------------------------------------

%% 4a | Read in lisst.ini variables
cfg.inst.ini = ini2struct(fullfile(cfg.path.dir_inst,'lisst.ini'));

%% 4b | Read instrument specifics
% read InstrumentData.txt - file from Sequoia with info based on factory
% configuration
cfg.inst.InstrumentData = readtable(fullfile(cfg.path.dir_inst,'InstrumentData.txt'),'Delimiter',',','FileType','text');
cfg.inst.VCC     = cfg.inst.InstrumentData.Var4;
cfg.inst.LISSTsn = cfg.inst.InstrumentData.Var1;
cfg.inst.typeX   = char(cfg.inst.InstrumentData.Var2); 

%% 4c | Determine type for processing size ranges
% Set a catch in case user defined type and type read from
% InstrumentData.txt is different. For example, this could be the case if
% the instrument type was changed as was done in June 2018 for LISST-Deep
% sn4025. 
if isfield(cfg.inst,'type')
  userdefined_type = cfg.inst.type;
end
switch cfg.inst.typeX
  case {'a', 'A'}; cfg.inst.type = 1; % TYPE = 1 for type A (5-500 µm)
  case {'b', 'B'}; cfg.inst.type = 2; % TYPE = 2 for type B (1.25-250 µm size range)
  case {'c', 'C'}; cfg.inst.type = 3; % TYPE = 3 for type C (2.5-500 µm size range)
  case {'d', 'D'}; cfg.inst.type = 4; % TYPE = 4 for FLOC   (7.5-1500 µm size range)
end
if exist('userdefined_type','var') && ~isequal(userdefined_type,cfg.inst.type)
  fprintf('**CAUATION** Instrument type from file is different than user input!\n')
  fprintf('Stopping here in keyboard  mode\n')
  keyboard
end
%% 4d | determine if the data are to be inverted in LISST-100/LISST-100X format (32 size bins)
x = char(cfg.inst.InstrumentData.Var5);
if strcmp(x,'X')
  cfg.inst.X    = 1;  % X = 1 if data were obtained with LISST-100X
else
  cfg.inst.X    = 0;  % X = 0 if data were NOT obtained with LISST-100X
end

%% 4e | Add instrument type that corresponds to compute_mean.m definitions
% type should be defined in [cfg.project]_config.m script but is updated
% here based on the model defined in cfg.proc_options.rand. 
% 1 for type A (discontinued)        | TYPE = 1 for type A (5-500 µm)
% 2 for type B                       | TYPE = 2 for type B (1.25-250 µm size range)
% 3 for type C                       | TYPE = 3 for type C (2.5-500 µm size range)
% 4 for FLOC (discontinued)          | TYPE = 4 for FLOC   (7.5-1500 µm size range)
% 21 for randomly shaped type B
% 31 for randomly shaped type C.
if cfg.inst.type == 2 && cfg.proc_options.rand == 1
  cfg.inst.type2 = 21;
elseif cfg.inst.type == 3 && cfg.proc_options.rand == 1
  cfg.inst.type2 = 31;
else
  cfg.inst.type2 = cfg.inst.type;
end

%% 5  | Quality control 
%% 5a | Choose whether to run through qc tests automatically, or prompt user for each cast
cfg.qaqc_options.expert_qc = 0;
fprintf('  Run through manual expert QC?\n')
fprintf('   <0> ONLY auto QC     - no user prompting \n')
fprintf('   <1> ALSO expert QC - step through casts \n')
cfg.qaqc_options.expert_qc = input('   Enter choice: ');

%% 5b | Define flags
% Code | Value | Definition
%  1   | Good | Passed documented QC tests
%  2   | Not evaluated, not available, or unknown | Used for data when no QC test was performed, or the information on quality is not available
%  3   | Questionable | Failed non-critical documented metric or subjective test
%  4   | Bad | Failed critical documented QC test(s) or as assigned by the data provider
%  5   | Estimate | Cell value were interpolated, extrapolated, or otherwise estimated
%  6   | Below detection limit | Value is below the detection limit of the analytical methods applied
%  7   | Outside instrument specification | Value is within the extended measuring range
%  9   | Missing Data | Used as placeholder when data are missing
cfg.qaqc_options.quality_flag = struct();
cfg.qaqc_options.quality_flag.good                    = 1;
cfg.qaqc_options.quality_flag.not_evaluated           = 2;
cfg.qaqc_options.quality_flag.questionable            = 3;
cfg.qaqc_options.quality_flag.bad                     = 4;
cfg.qaqc_options.quality_flag.estimate                = 5;
cfg.qaqc_options.quality_flag.below_detection_limit   = 6;
cfg.qaqc_options.quality_flag.outside_instrument_spec = 7;
cfg.qaqc_options.quality_flag.missing_data            = 9;

%% 5c | Define automatic quality control tests
%  < LISST-100XUsersManualVersion500.pdf >
% cfg.qaqc_options.test is a structure where each field is a unique QC test that
% is run in LISST_data_qaqc_auto.m.
% For example: bad_laser_ref is a QC test with fields....
%  method    | executable QC test. I.e. wflagged = eval(bad_laser_ref.method);
%  reason    | description of QC test
%  flag_type | executable QC flag. I.e. data.quality_flag(wflagged) = cfg.qaqc_options.quality_flag.(bad_laser_ref.flag_type);

cfg.qaqc_options.test.bad_laser_ref.method    = 'data_proc.laserPower <= 0.02;'; 
cfg.qaqc_options.test.bad_laser_ref.reason    = 'Laser Reference should never be 0 (or very close to 0, e.g. 0.02mW). If this is the case the laser is most likely dead. In this case only your pressure data and temperature data will be any good. The laser must be replaced. Contact Sequoia or your local distributor for instructions on what to do.';
cfg.qaqc_options.test.bad_laser_ref.flag_type = 'bad'; 

cfg.qaqc_options.test.transmission_bad1.method    = 'data_proc.tau <= 0 | data_proc.tau >= 1';
cfg.qaqc_options.test.transmission_bad1.reason    = 'The transmission must be a number between 0 and 1. It is physically impossible for the transmission to be negative or larger than 1 (one). If transmission shows up as being  larger than 1 (one), then your measurement is most likely taken in very clear water and/or you have a bad zscat measurement obtained with dirty water and/or dirty windows.';
cfg.qaqc_options.test.transmission_bad1.flag_type = 'bad';

cfg.qaqc_options.test.transmission_bad2.method    = 'data_proc.tau <= 0.10;';
cfg.qaqc_options.test.transmission_bad2.reason    = 'If your transmission values are < 0.10 (or 10%), the water is too turbid. Disregard these data.';
cfg.qaqc_options.test.transmission_bad2.flag_type = 'bad';

cfg.qaqc_options.test.transmission_bad3.method    = 'data_proc.tau >= 0.995;';
cfg.qaqc_options.test.transmission_bad3.reason    = 'If your transmission values are > 0.995 (or 99.5%), the water is too clear. Disregard these data.';
cfg.qaqc_options.test.transmission_bad3.flag_type = 'bad';

cfg.qaqc_options.test.transmission_questionable1.method    = 'data_proc.tau >= 0.10 & data_proc.tau <= 0.30;';
cfg.qaqc_options.test.transmission_questionable1.reason    = 'Be wary of data collected at transmission values between 0.30 and 0.10 -- generally decreasing data quality as the transmission decreases below 0.30 (30%).';
cfg.qaqc_options.test.transmission_questionable1.flag_type = 'questionable';

cfg.qaqc_options.test.transmission_questionable2.method    = 'data_proc.tau >= 0.980 & data_proc.tau <= 0.995;';
cfg.qaqc_options.test.transmission_questionable2.reason    = 'Be wary of data collected at transmission values between 0.98 and 0.995 (or 98%-99.5%) due to low signal-to-noise ratio. This means the data may have a lot of noise in them, but can most likely still be used.';
cfg.qaqc_options.test.transmission_questionable2.flag_type = 'questionable';

% simple QC based on number distribution from InLineAnalysis/processLISST.m
cfg.qaqc_options.test.negative_number_distribution.method    = 'any(data_proc.PSD_DNSD < 0,2)';
cfg.qaqc_options.test.negative_number_distribution.reason    = 'negative particle size distribution';
cfg.qaqc_options.test.negative_number_distribution.flag_type = 'bad';

end %% MAIN FUNCTION