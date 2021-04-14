function [cfg,data_raw, meta_raw] = LISST_read_raw_data(cfg)
%FUNCTION LISST_read_raw_data
%
%  Syntax:
%    [cfg, data_raw, meta_raw] = LISST_read_raw_data(cfg)
%  
%  Description:
%    Finds cruise binary files based on cfg.path.dir_raw and then allocates
%    which background scatterfile to use in getscat.m function. Background
%    scatterfile can be a single file (i.e. an average of all zscat over
%    the cruise), or can be individual zscat files. 
%
%  Refereces:
%    http://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/
%
%  Notes:
%    To do: review and incorporate http://www.sequoiasci.com/article/how-to-compute-the-mean-particle-diameter-from-a-lisst-volume-distribution-2/
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
dbstop if error
%% 0 | Set up base structure

if exist(cfg.path.file_raw,'file')
  fprintf('Loading %s\n',cfg.path.file_raw);
  load(cfg.path.file_raw,'data_raw','meta_raw','cfg')
  return
end

%% 1 | Initialize data and meta structures
data_raw = struct(); % structure containing data
meta_raw = struct(); % structure containing the same fields as data with variable name descriptions and units

%% 2 | Find data_raw filename
% Define basic structure that will contain file data
rawfiles = dir(fullfile(cfg.path.dir_raw,'*.DAT'));
rawfiles(contains({rawfiles.name},'._')) = [];
% remove empty files
bytes = [rawfiles.bytes];
idx_zero = find(bytes == 0);
rawfiles(idx_zero) = [];

% store number of files
num_rawfiles = numel(rawfiles);
if cfg.testing && num_rawfiles > 1
  num_rawfiles = 10;
end
% populate filenames
cfg.proc_options.datfile = {};
for nf = 1:num_rawfiles
  cfg.proc_options.datfile{nf} = rawfiles(nf).name;   % name of data_raw file
  cfg.proc_options.zscfile{nf} = '';                  % default to blank
end

%% 3 | Read raw data files (binary data_raw files)
fprintf('Reading binary data_raw data\n')
% raw (binary/data_raw) files are associated background (zscat) files 
% populate filenames
% column names of dat files (from Sequoia manual)
vars_1to32  = strcat({'lightintensity_detector'}, num2str([1:32]'));
vars_33to40 = {'laser_transmission' 'voltage_rawcounts' 'ext_aux' 'laser_reference' 'pressure_rawcounts' 'temperature_100th_degrees' 'day*100+hour' 'minutes*100+seconds'}';

for nraw = 1:num_rawfiles
  fprintf('Reading %s\n',cfg.proc_options.datfile{nraw});
  data = tt2mat(fullfile(cfg.path.dir_raw,cfg.proc_options.datfile{nraw}),40);  
  if nraw == 1
    data_raw.data    = data;
    data_raw.datfile = repmat(cfg.proc_options.datfile(nraw),size(data,1),1);
    data_raw.zscfile = repmat(cfg.proc_options.zscfile(nraw),size(data,1),1);
  else
    data_raw.data    = [data_raw.data; data];
    data_raw.datfile = [data_raw.datfile; repmat(cfg.proc_options.datfile(nraw),size(data,1),1)];
    data_raw.zscfile = [data_raw.zscfile; repmat(cfg.proc_options.zscfile(nraw),size(data,1),1)];
  end
end

% update meta
meta_raw.datfile.name = 'binary .DAT file offloaded from LISST instrument';
meta_raw.zscfile.name = 'zscat (background) file, typically obtained using the Windows SOP';
meta_raw.data.name    = 'raw data file, converted from binary';
meta_raw.data.fields  = [erase(vars_1to32,' '); vars_33to40];

%% 4 | Remove any duplicate entries and sort by date
% first calculate datenum
date = datenumfromdata(cfg.year,data_raw.data(:,39:40));
[~,uidx] = unique(date(isfinite(date))); % This will also sort by date
fields = fieldnames(data_raw);
for nf = 1:numel(fields)
  data_raw.(fields{nf}) = data_raw.(fields{nf})(uidx,:);
end

%% 5 | Save data to .mat file
if ~cfg.testing % only save if not in testing mode
  fprintf('Saving processed LISST data to %s\n',cfg.path.file_raw)
  save(cfg.path.file_raw,'cfg','data_raw','meta_raw','-v7.3');
end


end %% MAIN FUNCTION 