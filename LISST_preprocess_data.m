function [cfg, data_pre, meta_pre] = LISST_preprocess_data(cfg, data_raw, meta_raw)
%% function LISST_processing_workflow
%  Syntax: 
%     [cfg, data_pre, meta_raw] = LISST_preprocess_data(cfg, data_raw, meta_raw)
%
%  Description:
%    "Preprocesses" LISST data. 
%    1. Defines new structures "meta_pre" and "data_pre"
%    2. Calculates date (matlab's datenum format) 
%    3. Calculates temperatuer and depth based on lisst.ini file contents
%    4. Reads in CTD data ** subsequent steps much easier if this is upcast
%       and downcast data with true timestamp!
%    5. Matches LISST profile to CTD cast #
%    6. Corrects LISST depth and time based on CTD data
%    7. Identifies limits of downcast (decending) in LISST profile
%    8. Saves data to cfg.path.file_pre
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

%% 1 | Initialize data and meta structures
data_pre = data_raw; % structure containing data
meta_pre = meta_raw; % structure containing the same fields as data with variable name descriptions and units

%% 2 | Calculate date
fprintf('Calculating date\n')
data_pre.date     = datenumfromdata(cfg.year,data_pre.data(:,39:40));
if strcmp(cfg.project,'LISST_sn4041_2018_ASGARD_SKQ201813S')
  idx_2019 = year(data_pre.date) == 2019;
  data_pre.date(idx_2019) = data_pre.date(idx_2019) - datenum(1,0,0);
end
data_pre.datetime = datetime(data_pre.date,'ConvertFrom','datenum');

% update meta
meta_pre.date.name     = 'MATLAB datenum';
meta_pre.date.unit     = 'Number of days since 0-Jan-0000';
meta_pre.datetime.name = 'MATLAB datetime';
meta_pre.datetime.unit = 'scalar datetime array corresponding to the date and time';

%% 3 | Calculate temperature and depth
fprintf('Calculating temperature and depth\n')
% First, load necessary offsets and scales from ini file
instfield   = ['instrument' num2str(cfg.inst.LISSTsn)];
tempScale   = str2double(cfg.inst.ini.(instfield).hk5scale); % HK5Scale: temperature multiplier for calibration of temperature.
tempoffset  = str2double(cfg.inst.ini.(instfield).hk5off);   % HK5Off:   temperature pffset for calibration of temperature.
depthScale  = str2double(cfg.inst.ini.(instfield).hk4scale); % HK4Scale: pressure scales for calibration of pressure sensor. Look it up in LISST.INI for your serial #.
depthoffset = str2double(cfg.inst.ini.(instfield).hk4off);   % HK4Off:   pressure offsets for calibration of pressure sensor. Look it up in LISST.INI for your serial #.
% Calculate temperature and depth
data_pre.temp  = data_pre.data(:,38)*tempScale+tempoffset;
data_pre.depth = data_pre.data(:,37)*depthScale+depthoffset;
  
% update meta
meta_pre.temp.name  = 'temperature measured in endcap (therefore may have signficiant lag)';
meta_pre.temp.unit  = 'degC';
meta_pre.depth.name = 'depth calibrated using factory supplied constants from LISST.INI file';
meta_pre.depth.unit = 'm';


%% 4 | Read CTD data
if exist(cfg.path.file_ctddata,'file')
  fprintf('Loading CTD data from %s\n',cfg.path.file_ctddata)
  load(cfg.path.file_ctddata)
elseif exist(fullfile(cfg.path.dir_ctd,[cfg.project '_CTD.mat']),'file')
  cfg.path.file_ctddata = fullfile(cfg.path.dir_ctd,[cfg.project '_CTD.mat']);
  load(cfg.path.file_ctddata)
else
  ctd = read_ctd_data_by_type(cfg.ctd_type,cfg.path.dir_ctd);
end
if ~isfield(ctd,'station')
  fprintf('ADD STATION!\n')
  keyboard
end
% Calculate depth from pressure and latitude
if ~isfield(ctd,'depth')
  ctd.depth = -1*gsw_z_from_p(ctd.press,ctd.lat);
end

%% 5 | Match LISST profile to CTD sequential cast number
[cfg, data_pre, meta_pre] = LISST_match_ctd_cast(cfg, data_pre, meta_pre, ctd);

%% 6 | Remove all casts with bad data
if strcmp(cfg.project,'LISST_sn4025_2019_NGA_TGX201909') || strcmp(cfg.project,'LISST_sn4025_2017_ASGARD_SKQ201709S')
  % Depth was not zero'd before deployment which resulted in a ~340m depth
  % offset, so any data below 340 m is not recoverable!
  ucasts = unique(data_pre.cast);
  bad_casts  = [];
  good_casts = [];
  for ncast = 1:numel(ucasts)
    idx = data_pre.cast == ucasts(ncast);
    if all(data_pre.depth(idx) == 0 | isnan(data_pre.depth(idx)))
      bad_casts  = [bad_casts;ucasts(ncast)];
    else
      good_casts = [good_casts;ucasts(ncast)];
    end
  end
  % if any casts were identified as bad... remove those casts
  if ~isempty(bad_casts)
    % Move raw files to new folder
    idx_rm_casts = ismember(data_pre.cast,bad_casts);
    move_datfiles = unique(data_pre.datfile(idx_rm_casts));
    zero_dir = fullfile(cfg.path.dir_raw,'_depth_all_zero');
    if ~isdir(zero_dir)
      mkdir(zero_dir)
    end
    
    for nrm = 1:numel(move_datfiles)
      if exist(fullfile(cfg.path.dir_raw,move_datfiles{nrm}),'file')
        movefile(fullfile(cfg.path.dir_raw,move_datfiles{nrm}),fullfile(zero_dir,move_datfiles{nrm}));
      elseif exist(fullfile(zero_dir,move_datfiles{nrm}),'file')
        fprintf('%s already in "_depth_all_zero" folder\n',move_datfiles{nrm})
      else
        fprintf('could not locate %s...\n',move_datfiles{nrm})
        keyboard
      end
    end
    
    % Remove bad casts
    idx_remove = ismember(data_pre.cast,bad_casts);
    t = struct2table(data_pre);
    t(idx_remove,:) = [];
    % Set data to NaN where depth == 0
    idx_remove = t.depth == 0;
    t(idx_remove,:) = [];
    
    data_pre = table2struct(t,'ToScalar',true);
    idx_remove = ismember(cfg.proc_options.cast,bad_casts);
    cfg.proc_options.cast(idx_remove) = [];
    cfg.proc_options.datfile(idx_remove) = [];
    cfg.proc_options.zscfile(idx_remove) = [];
    
  else
    % Set data to NaN where depth == 0
    t = struct2table(data_pre);
    idx_remove = t.depth == 0;
    t(idx_remove,:) = [];
    data_pre = table2struct(t,'ToScalar',true);
  end
end

  
%% 6 | Correct time and depth lag
[cfg, data_pre, meta_pre] = LISST_correct_time_depth_lag(cfg,data_pre, meta_pre,ctd);

%% 7 | Limit data to downcast
[cfg, data_pre, meta_pre] = LISST_identify_downcast(cfg,data_pre,meta_pre,ctd);
keyboard
%% 8 | Save data_pre and meta_pre structures 
fprintf('Saving preprocessed data to file: %s\n',cfg.path.file_pre)
save(cfg.path.file_pre,'cfg','data_pre','meta_pre');
end %% MAIN FUNCTION