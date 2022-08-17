function [cfg, data_proc, meta_proc] = LISST_process_data(cfg,data_pre, meta_pre)
%FUNCTION LISST_process_data
%
%  Syntax:
%    [cfg, data_proc] = LISST_process_data(cfg,data_proc)
%
%  Description:
%    Processes LISST data files and stores entire cruise dataset to
%    data_proc structure. Variable description and units stored in meta_proc
%    structure.
%
%  Refereces:
%    http://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/
%    LISST-Deep-Users-Manual-May-2013.pdf
%
%  Notes:
%    All units with micro prefix "µ" noted as "u"
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
dbstop if error
close all
%% 0 | Set up basic information
% save with information about what background file and inversion were used
if cfg.proc_options.rand == 1 % non-spherical inversion model
  model_str = 'nonspherical';
elseif cfg.proc_options.rand == 0 % spherical inversion model
  model_str = 'spherical';
else
  fprintf('Improper inversion model selected, can only be 1 (non-spherical) or 0 (spherical) currently\n')
  keyboard
end
cfg.path.file_proc = fullfile(cfg.path.dir_proc,[cfg.project '_proc_' model_str '.mat']);  

%% 1 | Initialize data strcutures
data_proc = data_pre;
meta_proc = meta_pre;

%% 2 | Limit to downcasts only
idx_downcasts = data_proc.profile_index == 1;
fields = fieldnames(data_proc);
for nf = 1:numel(fields)
  data_proc.(fields{nf}) = data_proc.(fields{nf})(idx_downcasts,:);
end
  
%% 3 | Select which background to use
if isfield(cfg,'zscfile') && numel(cfg.zscfile) == 1 && exist(fullfile(cfg.path.dir_zsc,char(cfg.zscfile)),'file')
  data_proc.zscfile(:)        = cfg.zscfile;
  cfg.proc_options.zscfile(:) = cfg.zscfile;
elseif numel(unique(cfg.proc_options.zscfile)) == 1 && ~isempty(char(cfg.proc_options.zscfile(1))) && exist(fullfile(cfg.path.dir_zsc,char(cfg.proc_options.zscfile(1))),'file')
  data_proc.zscfile(:)        = cfg.zscfile;
  cfg.proc_options.zscfile(:) = cfg.zscfile;
else
  [cfg, data_proc] = LISST_select_background_scatterfiles(cfg,data_proc);
end

keyboard
%% 4 | Pull out size bin information for LISST instrument   
% Depends on instrument, inversion type and laser color
cfg.inst.bins = LISST_bin_sizes(cfg.inst.type2);
num_szbin = numel(cfg.inst.bins.bin_size);
num_scans = size(data_proc.data,1);
% Create string
bins.name   = 'bin size based on instrument type and inversion model';
bins.unit   = 'um';
bins.nlimits = [];
bins.slimits = {};
for nb = 1:numel(cfg.inst.bins.low_limit)
  bins.nlimits(nb,1:2) = [cfg.inst.bins.low_limit(nb) cfg.inst.bins.upp_limit(nb)];
  bins.slimits(nb,1)     = {[num2str(cfg.inst.bins.low_limit(nb)) '-' num2str(cfg.inst.bins.upp_limit(nb)) 'um']};
end

%% 5 | Pull out angles for LISST 100X VSF measurements in water
% Read in angles (in degrees) for VSF measurements
% <https://www.sequoiasci.com/article/angles-for-lisst-100x-vsf-measurement-in-water/>
angles = readtable(fullfile(cfg.path.toolbox,'lisst_utility','angles_for_vsf_measurement_in_water.csv'));
angles = flipud(angles); % change so SizebinNumber goes from 1:32 rather than RingNmber
% Only store the angles for the specific type of lisst instrument
gd = contains(angles.Properties.VariableNames,{'RingNumber' 'SizebinNumber' ['type' cfg.inst.typeX]});
angles = angles(:,gd);

%% 6 | Define constants
% needed for invert script
SHARPEN     = 0; % 1 causes the routine to check if the size distribution is narrow and, if so, increases the number of inversions. Use this setting if you expect a narrow size distribution (e.g. if you are analyzing narrow-size standard particles).
WAITBARSHOW = 0; % 1 if user wants a waitbar to show during processing in order to keep track of progress.
% needed for get_bf_v2 script
instname = ['instrument' num2str(cfg.inst.LISSTsn)];
HK3Scale    = str2double(cfg.inst.ini.(instname).hk3scale); % HK3Scale: mW laser power entering water per digital count in factory zscat variable 36; Look it up in LISST.INI for your serial #.
depthScale  = str2double(cfg.inst.ini.(instname).hk4scale); % HK4Scale: pressure scales for calibration of pressure sensor. Look it up in LISST.INI for your serial #.
depthoffset = str2double(cfg.inst.ini.(instname).hk4off);   % HK4Off: pressure offsets for calibration of pressure sensor. Look it up in LISST.INI for your serial #.
% needed for compute_mean script
name_cols = {'tot_vol_concentration' 'mean_size' 'std' 'tau' 'D10' 'D16' 'D50' 'D60' 'D84' 'D90' ...
  'hazen_u_coef' 'surf_area' 'silt_dens' 'silt_vol'};
unit_cols = {'uL/L' 'um' 'um' 'nodim' 'um' 'um' 'um' 'um' 'um' 'um' 'nodim' 'cm2/l' 'nodim' 'um'};
desc_cols = {'Total Volume Concentration' 'Mean Size' 'Standard Deviation' 'transmission' 'D10' 'D16' 'D50' 'D60' 'D84' 'D90' ...
  'D60/D10 (the Hazen uniformity coefficient)' ...
  'Surface area' 'Silt Density - The volume concentration ratio of silt particles to the total volume concentration'...
  'Silt Volume - the volume concentration of all particles < 64 um.'};
bin_size = cfg.inst.bins.bin_size;
if ~isequal(size(bin_size,2),num_szbin)
  bin_size = bin_size';
end
% needed for converting cscat to VSF in scientific units
cfg.inst.path_m = cfg.inst.path_length./100; % cm to meter
cfg.inst.phi    = 1/6;                       % Fraction of circle covered by detectors -- is this where the units of 1/m/sr comes from?
try   sq_size = cfg.inst.bins.ds(2:end).^2 - cfg.inst.bins.ds(1:end-1).^2;
catch sq_size = cfg.inst.bins.upp_limit(1:end).^2 - cfg.inst.bins.low_limit(1:end).^2; end
if ~isequal(size(sq_size,2),num_szbin); sq_size = sq_size'; end

%% 7 | Initialize variables and add metadata
% Update this as more variables are added!
data_proc.scat  = nan(num_scans,num_szbin);
data_proc.cscat = nan(num_scans,num_szbin);
data_proc.tau   = nan(size(data_proc.date));
meta_proc.scat.name  = 'raw scattering signature in digital counts (i.e. NOT corrected with the ringarea file!!)';
meta_proc.scat.unit  = 'digital counts';
meta_proc.tau.name   = 'optical transmission';                                   % http://www.sequoiasci.com/article/optical-transmission-what-is-it/
meta_proc.cscat.name = 'corrected scattering, corrected with the ringarea file';
meta_proc.cscat.unit = 'digital counts';

data_proc.uncal_vd  = nan(num_scans,num_szbin);
meta_proc.uncal_vd.name   = 'uncalibrated volume distribution';
meta_proc.uncal_vd.unit   = 'uL/L';
meta_proc.uncal_vd.method = 'Sequoias proprietary invert.p';

meta_proc.dias.name       = 'midpoint of size bins';
meta_proc.dias.unit       = 'um';
meta_proc.dias.method     = 'the midpoint of the size bins for the 8 / 32 size classes for  the appropriate instrument, inversion type and laser color calculated with Sequoias proprietary routine invert.p';

data_proc.VSD = nan(num_scans,num_szbin);
meta_proc.VSD.name   = 'Volume size distribution';
meta_proc.VSD.unit   = 'uL/L'; % 1uL = 1mm^3, so uL/L=mm^3/L=ppm
meta_proc.VSD.method = 'Volume distribution calibrated with the Volume Conversion Constant gives the volume concentration of particles in the different size classes';
meta_proc.VSD.unit2  = '1uL = 1mm^3, so uL/L == mm^3/L == ppm'; % add note about unit equality/conversion
meta_proc.VSD.bins   = bins;


data_proc.PSD  = nan(num_scans,num_szbin);
meta_proc.PSD.name  = 'particle size distribution';
meta_proc.PSD.unit  = '#/L';
meta_proc.PSD.bins  = bins;


data_proc.PSD_DNSD  = nan(num_scans,num_szbin);
meta_proc.PSD_DNSD.name  = 'Differential number size distribution';
meta_proc.PSD_DNSD.unit  = '#/m^3/um';
meta_proc.PSD_DNSD.bins  = bins;

data_proc.PSD_DVSD  = nan(num_scans,num_szbin);
meta_proc.PSD_DVSD.name  = 'Differential volume size distribution';
meta_proc.PSD_DVSD.unit  = 'uL/m^3/um';  % 1ppm = 1mm^3/L
meta_proc.PSD_DVSD.bins  = bins;

for nc = 1:numel(name_cols)
  if nc == 4; continue; end % tau already a variable
  data_proc.(name_cols{nc}) = nan(size(data_proc.date));
  meta_proc.(name_cols{nc}).unit = unit_cols{nc};
  meta_proc.(name_cols{nc}).name = desc_cols{nc};
end

data_proc.bf               = nan(size(data_proc.date));
data_proc.beam_attenuation = nan(size(data_proc.date));
data_proc.absorption       = nan(size(data_proc.date));
data_proc.laserPower       = nan(size(data_proc.date));
meta_proc.bf.name               = 'forward scattering coefficient up to max range covered by LISST';
meta_proc.bf.unit               = '1/m';
meta_proc.beam_attenuation.name = 'beam attenuation, calculated using path length';
meta_proc.beam_attenuation.unit = '1/m';
meta_proc.absorption.name       = 'absorption = beam attentuation - forward scattering coefficient';
meta_proc.absorption.unit       = '';
meta_proc.laserPower.name       = 'laser power entering water, in units of mW';
meta_proc.laserPower.unit       = 'mW';

data_proc.VSF  = nan(num_scans,num_szbin);
meta_proc.VSF.name = 'Volume Scattering Function';
meta_proc.VSF.unit = '1/m/sr';
meta_proc.VSF.angles = angles;
meta_proc.VSF.bins   = bins;

data_proc.nVSF = nan(num_scans,num_szbin);
meta_proc.nVSF.name   = 'Normalized Volume Scattering Function, VSF divided by beam attenutation';
meta_proc.nVSF.unit   = '1/sr';
meta_proc.nVSF.angles = angles;
meta_proc.nVSF.bins   = bins;

%% 8 | Add header information for the differential fields
% These are required for SeaBASS, so define here. 
%/PSD_bin_size_method=specify method used to select the nominal bin size such as arithmetic_mean or?geometric_mean
%/PSD_bin_size_boundaries=please provide a?comma-separated?list with the bin size boundaries in increasing order, e.g., 5,9.5,15,20 ...
cfg.header.PSD_bin_size_boundaries = erase(strjoin(cellstr(num2str(unique([cfg.inst.bins.low_limit;cfg.inst.bins.upp_limit])))',','),' ');
switch cfg.inst.type2
  case {1 2 3 4} % LISST type  A, B, C, and FLOC
    % mid_point = sqrt(lower_limit.*upper_limit);
    cfg.header.PSD_bin_size_method = 'square root of the product of the lower and upper bin limits';
  case {21 31} %
    % Defined in LISST_bin_sizes.m, adapted Sequoia's compute_mean.m
    % <https://www.sequoiasci.com/article/compute_mean-m/>
    cfg.header.PSD_bin_size_method = 'Defined by Sequoia compute_mean.m';
end

%% 9 | Range
% Calculate parameters using the highest and lowest backgrounds, the range
% of those calculations can give an estimate of uncertainty
if cfg.proc_options.calculate_range && isfield(cfg,'range')
  range_variables = {'scat' 'cscat' 'tau' 'uncal_vd' 'VSD' 'tot_vol_concentration' 'mean_size' 'PSD' 'PSD_DNSD' 'PSD_DVSD'}; % Include uncertainty estimation for these parameters
  for nv = 1:numel(range_variables)
    var = range_variables{nv};
    var_rng1 = [var '_range1']; % Lowest background (that is NOT factory) and highest background
    var_rng2 = [var '_range2']; % factory background and highest background - this will be higher (unless something really unexpected happened)
    
    data_proc.(var_rng1) = nan(size(data_proc.(var)));
    meta_proc.(var_rng1) = meta_proc.(var);
    meta_proc.(var_rng1).name   = ['range of ' meta_proc.(var_rng1).name];
    meta_proc.(var_rng1).zscat  = cfg.range.range1.zscfile;
    meta_proc.(var_rng1).source = ['absolute difference between ' var ' dervied using minimum background from ' cfg.range.range1.zscfile{1} ' and maximum background from ' cfg.range.range1.zscfile{2}];
    
    data_proc.(var_rng2) = nan(size(data_proc.(var)));
    meta_proc.(var_rng2) = meta_proc.(var);
    meta_proc.(var_rng2).name   = ['range of ' meta_proc.(var_rng2).name];
    meta_proc.(var_rng2).zscat  = cfg.range.range2.zscfile;
    meta_proc.(var_rng2).source = ['absolute difference between ' var ' dervied using factory background from ' cfg.range.range2.zscfile{1} ' and maximum background from ' cfg.range.range2.zscfile{2}];
  end
end

%% ************************************************************************
%% 10 | LOOP THROUGH UNIQUE BACKGROUNDS AND PROCESS RAW DATA
%% ************************************************************************
% Although all the data is in array format, still need to loop through
% unique background files (sometimes correlates to casts) to process.
% This way of processing will hopefully allow for processing of moored
% sensor data as well as data from sensors mounted on a rossette
background_files = unique(data_proc.zscfile,'stable');
for nzscat = 1:numel(background_files)
  idx_zscat = find(strcmp(data_proc.zscfile, background_files{nzscat}));
  fprintf('Processing data with background #%d of %d: %s....................................................\n',nzscat,numel(background_files),background_files{nzscat})
  
  %% I    | Process raw data to get scattering and transmission
  % The first step is to convert the binary .DAT file with scattering data into
  % corrected scattering data (cscat) using the background (zscat) file and
  % ring area file.
  fprintf(' Converting raw data into corrected scattering data\n')
  try
    zscfile = fullfile(cfg.path.dir_zsc,cfg.proc_options.zscfile{nzscat});
    [data_proc.scat(idx_zscat,:), data_proc.tau(idx_zscat), ~, ~, data_proc.cscat(idx_zscat,:)] = getscat_v2(data_proc.data(idx_zscat,:),zscfile,cfg.inst.X,cfg.path.file_ringarea);
    % % commented out April 2021 - not necessary
    % data_proc.zsc(idx_zscat,:)  = repmat(zsc,numel(idx_zscat),1); 
  catch
    fprintf('  problem with getscat with background file: %s\n',zscfile)
    keyboard
  end
  if cfg.proc_options.calculate_range
    try
      range_zscfile.fac = fullfile(cfg.path.dir_zsc,cfg.range.range2.zscfile{1}); % factory background
      range_zscfile.min = fullfile(cfg.path.dir_zsc,cfg.range.range1.zscfile{1}); % lowest measured background, other than factory
      range_zscfile.max = fullfile(cfg.path.dir_zsc,cfg.range.range1.zscfile{2}); % highest measured background
      [scat.fac, tau.fac, ~, ~, cscat.fac] = getscat_v2(data_proc.data(idx_zscat,:),range_zscfile.fac,cfg.inst.X,cfg.path.file_ringarea);
      [scat.min, tau.min, ~, ~, cscat.min] = getscat_v2(data_proc.data(idx_zscat,:),range_zscfile.min,cfg.inst.X,cfg.path.file_ringarea);
      [scat.max, tau.max, ~, ~, cscat.max] = getscat_v2(data_proc.data(idx_zscat,:),range_zscfile.max,cfg.inst.X,cfg.path.file_ringarea);
      data_proc.scat_range1(idx_zscat,:)  = abs(scat.max  - scat.min);
      data_proc.scat_range2(idx_zscat,:)  = abs(scat.max  - scat.fac);
      data_proc.cscat_range1(idx_zscat,:) = abs(cscat.max - cscat.min);
      data_proc.cscat_range2(idx_zscat,:) = abs(cscat.max - cscat.fac);
      data_proc.tau_range1(idx_zscat)     = abs(tau.max   - tau.min);
      data_proc.tau_range2(idx_zscat)     = abs(tau.max   - tau.fac);
      % idx = data_proc.cast == 30; % choose arbritrary cast to check
      % figure; subplot(2,1,1); errorbar(data_proc.tau(1:1000:end),data_proc.tau_unc(1:1000:end))
      % subplot(2,1,2); plot(data_proc.tau(idx),data_proc.depth(idx),'k*');
      % hold on;        plot(data_proc.tau_unc(idx),data_proc.depth(idx),'r.');
      % figure; errorbar(data_proc.tau(idx),data_proc.depth(idx),data_proc.tau_unc(idx));
    catch
      fprintf('  problem calculating uncertainty')
      fprintf('  background file %s or %s or %s\n',range_zscfile.fac, range_zscfile.min, range_zscfile.max)
      cfg.proc_options.calculate_range = 0;
    end
  end
  
  
  %% II   | Get uncalibrated volume distribution and midpoint of size bins in microns
  fprintf(' Calculating uncalibrated volume distribution and midpoint of size bins\n')
  [data_proc.uncal_vd(idx_zscat,:), dias] = invert(data_proc.cscat(idx_zscat,:), cfg.inst.type, cfg.inst.ST, cfg.proc_options.rand, SHARPEN, cfg.inst.green, WAITBARSHOW);
  
  % Define dias in cfg structure and check that it is as expected. 
  % dias = midpoint of the size bins for the 8 or 32 size classes for the
  % appropriate instrument, inversion type and laser color
  if nzscat == 1
    cfg.proc_options.dias = dias;
    if ~isequal(num2str(cfg.inst.bins.mid_point','%.4f'),num2str(dias,'%.4f'))
      % mid size of bins are not the same as those used in invert function??
      for ns = 1:num_szbin
        fprintf('%f  %f  %f\n',cfg.inst.bins.mid_point(ns),cfg.proc_options.dias(ns),cfg.inst.bins.mid_point(ns)-dias(ns))
      end
      fprintf('ERROR** mid size of bin not same as expected??\n')
      keyboard
    end
  end
  % Calculate associated uncertainty 
  if cfg.proc_options.calculate_range && ismember('uncal_vd',range_variables)
    % Use the corrected scatter derived from minimum and maximum
    % backgrounds, then take the difference to get an uncertainty estimate.
    [uncal_vd.fac, ~] = invert(cscat.fac, cfg.inst.type, cfg.inst.ST, cfg.proc_options.rand, SHARPEN, cfg.inst.green, WAITBARSHOW);
    [uncal_vd.min, ~] = invert(cscat.min, cfg.inst.type, cfg.inst.ST, cfg.proc_options.rand, SHARPEN, cfg.inst.green, WAITBARSHOW);
    [uncal_vd.max, ~] = invert(cscat.max, cfg.inst.type, cfg.inst.ST, cfg.proc_options.rand, SHARPEN, cfg.inst.green, WAITBARSHOW);
    data_proc.uncal_vd_range1(idx_zscat,:) = abs(uncal_vd.max - uncal_vd.min);
    data_proc.uncal_vd_range2(idx_zscat,:) = abs(uncal_vd.max - uncal_vd.fac);
  end
  
  %% III  | Get calibrated volume distribution, or VSD
  % next step is to invert the data to generate size distribution, which
  % we also call volume distribution since it is actually the volume
  % concentration of particles in the different size classes.
  % Convert the volume distribution into calibrated units, you must
  % divide uncal_vd by the Volume Conversion Constant (VCC)
  fprintf(' Calculating calibrated volume distribution\n')
  laser_ref = data_proc.data(idx_zscat,36);  % laster reference value during measurement
  if ~isfield(cfg.inst,'flref') % flref is the factory laser reference value for this instrument
    factory_zscat = load(cfg.path.file_factoryz);
    cfg.inst.flref = factory_zscat(36);
  end
  data_proc.VSD(idx_zscat,:) = vdcorr(data_proc.uncal_vd(idx_zscat,:), cfg.inst.VCC, cfg.inst.flref, laser_ref);
  % Calculate associated uncertainty
  if cfg.proc_options.calculate_range && ismember('VSD',range_variables)
    VSD.fac = vdcorr(uncal_vd.fac, cfg.inst.VCC, cfg.inst.flref, laser_ref);
    VSD.min = vdcorr(uncal_vd.min, cfg.inst.VCC, cfg.inst.flref, laser_ref);
    VSD.max = vdcorr(uncal_vd.max, cfg.inst.VCC, cfg.inst.flref, laser_ref);
    data_proc.VSD_range1(idx_zscat,:) = abs(VSD.max - VSD.min);
    data_proc.VSD_range2(idx_zscat,:) = abs(VSD.max - VSD.fac);
  end
  
  %% IV   | Calculate number size distribution, or PSD
  % Formulas adapted from Emmanuel Boss's processLISST.m data (InLineAnalysis package)
  fprintf(' Calculating number size distribution\n')
  data_proc.PSD(idx_zscat,:) = ( data_proc.VSD(idx_zscat,:) ./ (pi*cfg.proc_options.dias.^3/6) ).*10^9; % units: [uL/L] * [1/um^3] = [10^9um^3/L] * [1/um^3] = 10^9 [#/L]
  % Calculate ranges
  if cfg.proc_options.calculate_range && ismember('PSD',range_variables)
    PSD.fac = ( VSD.fac ./ (pi*cfg.proc_options.dias.^3/6) ).*10^9; % [#/L]
    PSD.min = ( VSD.min ./ (pi*cfg.proc_options.dias.^3/6) ).*10^9; % [#/L]
    PSD.max = ( VSD.max ./ (pi*cfg.proc_options.dias.^3/6) ).*10^9; % [#/L]
    data_proc.PSD_range1(idx_zscat,:) = abs(PSD.max - PSD.min);
    data_proc.PSD_range2(idx_zscat,:) = abs(PSD.max - PSD.fac);
  end
  
  %% V    | Calculate differential PSD & differential VSD
  % Resize the width of the size bins to include N rows
  bs = bin_size .* ones(numel(idx_zscat),1); % e.g. [1x32] to [Nx32]
  % Calculate differential fields by dividing by the width of the size bin.
  % Then convert per liter to per cubic meter 
  % [uL/L/um] --> [uL/m^3/um]  [#/L/um] *[1000L/m^3] = [#/m^3/um]
   % [#/L/um] --> [#/m^3/um].  [uL/L/um]*[1000L/m^3] = [uL/m^3/um]
  data_proc.PSD_DNSD(idx_zscat,:) = (data_proc.PSD(idx_zscat,:) ./ bs).*1000; % units: [#/L] / [um] *[1000L/m^3] = [#/m^3/um]  or [number/m^3/um]
  data_proc.PSD_DVSD(idx_zscat,:) = (data_proc.VSD(idx_zscat,:) ./ bs).*1000; % units: [uL/L]/ [um] *[1000L/m^3] = [uL/m^3/um] or [mm^3/m^3/um] 
  
  % Calculate ranges
  if cfg.proc_options.calculate_range && ismember('PSD_DNSD',range_variables)
    PSD_DVSD.fac = ( VSD.fac ./ bs).*1000;
    PSD_DVSD.min = ( VSD.min ./ bs).*1000;
    PSD_DVSD.max = ( VSD.max ./ bs).*1000;
    
    PSD_DNSD.fac = ( PSD.fac ./ bs).*1000;
    PSD_DNSD.min = ( PSD.min ./ bs).*1000;
    PSD_DNSD.max = ( PSD.max ./ bs).*1000;
    
    data_proc.PSD_DVSD_range1(idx_zscat,:) = abs(PSD_DVSD.max - PSD_DVSD.min); % range1 = max-min
    data_proc.PSD_DVSD_range2(idx_zscat,:) = abs(PSD_DVSD.max - PSD_DVSD.fac); % range2 = max-factory
    data_proc.PSD_DNSD_range1(idx_zscat,:) = abs(PSD_DNSD.max - PSD_DNSD.min); % range1 = max-min
    data_proc.PSD_DNSD_range2(idx_zscat,:) = abs(PSD_DNSD.max - PSD_DNSD.fac); % range2 = max-factory
  end
  
  %% VI   | Calculate mean particle size
  % variable and variable units in compute_mean.m notes
  fprintf(' Calculating mean particle size\n')
  variables_out = compute_mean(data_proc.VSD(idx_zscat,:), cfg.inst.type2, data_proc.tau(idx_zscat));
  % create new variables and fill with data from variables_out matrix
  for nc = 1:numel(name_cols)
    if nc == 4;continue;end % ignore because either NaN or tau (transmission) if tau was an inputend
    data_proc.(name_cols{nc})(idx_zscat) = variables_out(:,nc);
  end
  
  % Calculate ranges
  if cfg.proc_options.calculate_range && all(ismember({'tau' 'VSD'},range_variables))
    variables.fac = compute_mean(VSD.fac, cfg.inst.type2, tau.fac);
    variables.min = compute_mean(VSD.min, cfg.inst.type2, tau.min);
    variables.max = compute_mean(VSD.max, cfg.inst.type2, tau.max);
    % create new variables and fill with data from variables_out matrix
    for nc = 1:numel(name_cols)
      if nc == 4;continue;end % ignore because either NaN or tau (transmission) if tau was an input
      if ismember(name_cols{nc},range_variables)
        data_proc.([name_cols{nc} '_range1'])(idx_zscat) = abs(variables.max(:,nc) - variables.min(:,nc));
        data_proc.([name_cols{nc} '_range2'])(idx_zscat) = abs(variables.max(:,nc) - variables.fac(:,nc));
      end
    end
    clearvars variables.fac variables.min variables.max
  end
  
  %% VII  | Calculate forward scattering coefficient
  % This MATLAB script (http://www.sequoiasci.com/article/get_bf-m/) allows
  % the LISST-100 or LISST-100X user to compute the part of the forward
  % scattering coefficient, bf, covered by the angle ranges of the
  % LISST-100X. See the application note ‘Angles for LISST-100X VSF
  % measurement in WATER’ for details about the angles covered by a type B
  % and a type C LISST.
  fprintf(' Calculating forward scattering coefficient, beam attenuation and absorption\n')
  [data_proc.bf(idx_zscat), data_proc.beam_attenuation(idx_zscat), data_proc.absorption(idx_zscat), data_proc.laserPower(idx_zscat)] = get_bf_v2(data_proc.tau(idx_zscat), data_proc.data(idx_zscat,:), data_proc.cscat(idx_zscat,:), cfg.path.file_ringarea, HK3Scale, depthScale, depthoffset, cfg.inst.path_length, cfg.inst.X);
  
  % Calculate ranges
  if cfg.proc_options.calculate_range && any(ismember({'bf' 'beam_attenuation' 'absorption' 'laserPower'},range_variables))
    [bf.fac, beam_attenuation.fac, absorption.fac, laserPower.fac] = get_bf_v2(tau.fac, data_proc.data(idx_zscat,:), cscat.fac, cfg.path.file_ringarea, HK3Scale, depthScale, depthoffset, cfg.inst.path_length, cfg.inst.X);
    [bf.min, beam_attenuation.min, absorption.min, laserPower.min] = get_bf_v2(tau.min, data_proc.data(idx_zscat,:), cscat.min, cfg.path.file_ringarea, HK3Scale, depthScale, depthoffset, cfg.inst.path_length, cfg.inst.X);
    [bf.max, beam_attenuation.max, absorption.max, laserPower.max] = get_bf_v2(tau.max, data_proc.data(idx_zscat,:), cscat.max, cfg.path.file_ringarea, HK3Scale, depthScale, depthoffset, cfg.inst.path_length, cfg.inst.X);
    if ismember('bf',range_variables)
      data_proc.bf_range1(idx_zscat) = abs(bf.max - bf.min); % range2 = max-min
      data_proc.bf_range2(idx_zscat) = abs(bf.max - bf.fac); % range2 = max-factory
    end
    if ismember('beam_attenuation',range_variables)
      data_proc.beam_attenuation_range1(idx_zscat) = abs(beam_attenuation.max - beam_attenuation.min); % range2 = max-min
      data_proc.beam_attenuation_range2(idx_zscat) = abs(beam_attenuation.max - beam_attenuation.fac); % range2 = max-factory
    end
    if ismember('absorption',range_variables)
      data_proc.absorption_range1(idx_zscat) = abs(absorption.max - absorption.min); % range2 = max-min
      data_proc.absorption_range2(idx_zscat) = abs(absorption.max - absorption.fac); % range2 = max-factory
    end
    if ismember('laserPower',range_variables)
      data_proc.laserPower_range1(idx_zscat) = abs(laserPower.max - laserPower.min); % range2 = max-min
      data_proc.laserPower_range2(idx_zscat) = abs(laserPower.max - laserPower.fac); % range2 = max-factory
    end
    clearvars bf absorption laserPower
  end
  
  %% VIII | Convert VSF from counts (scat) to scientific units (1/m) and calculate scatter phase function (or normalized volume scattering function VSF/c)
  %http://www.sequoiasci.com/article/optical-volume-scattering-function-vsf-measurement-with-lissts/
  % Formulas adapted from Emmanuel Boss's processLISST.m data (InLineAnalysis package)
  fprintf(' Calculating Volume Scattering Function\n')
  data_proc.VSF(idx_zscat,:)  = data_proc.cscat(idx_zscat,:) ./ (ones(numel(idx_zscat),1) .* (pi * cfg.inst.phi * cfg.inst.path_m  * sq_size));
  data_proc.nVSF(idx_zscat,:) = data_proc.VSF(idx_zscat,:)   ./ data_proc.beam_attenuation(idx_zscat);
  
  % Calculate ranges
  if cfg.proc_options.calculate_range 
    if ismember('VSF',range_variables)
      VSF.fac  = cscat.fac ./ (ones(numel(idx_zscat),1) .* (pi * cfg.inst.phi * cfg.inst.path_m  * sq_size));
      VSF.min  = cscat.fac ./ (ones(numel(idx_zscat),1) .* (pi * cfg.inst.phi * cfg.inst.path_m  * sq_size));
      VSF.max  = cscat.fac ./ (ones(numel(idx_zscat),1) .* (pi * cfg.inst.phi * cfg.inst.path_m  * sq_size));
      data_proc.VSF_range1(idx_zscat) = abs(VSF.max - VSF.min); % range2 = max-min
      data_proc.VSF_range2(idx_zscat) = abs(VSF.max - VSF.fac); % range2 = max-factory
      if all(ismember({'nVSF' 'beam_attenuation'},range_variables))
        nVSF.fac = VSF.fac ./ beam_attenuation.fac;
        nVSF.min = VSF.min ./ beam_attenuation.min;
        nVSF.max = VSF.max ./ beam_attenuation.max;
        data_proc.nVSF_range1(idx_zscat) = abs(nVSF.max - nVSF.min); % range2 = max-min
        data_proc.nVSF_range2(idx_zscat) = abs(nVSF.max - nVSF.fac); % range2 = max-factory
        clearvars beam_attenuation nVSF VSF
      end % nVSF
    end % VSF
  end % ranges
end

%% 11 | Save data to .mat file
fprintf('Saving processed LISST data to %s\n',cfg.path.file_proc)
save(cfg.path.file_proc,'cfg','data_proc','meta_proc','-v7.3');

end %% MAIN FUNCTION