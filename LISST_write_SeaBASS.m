function LISST_write_SeaBASS(cfg,data_proc,meta_proc,data_grid)
%FUNCTION LISST_write_SeaBASS
%
%  Syntax:
%    LISST_write_SeaBASS(cfg,data_grid,meta)
%
%  Description:
%
%  Refereces:
%     https://seabass.gsfc.nasa.gov/wiki/Data_Submission
%     https://seabass.gsfc.nasa.gov/wiki/stdfields#Table%20of%20Field%20Names%20and%20Units
%
%  Notes:
%    Check file format with fcheck - will give exact details on any errors
%    or warnings generated. <https://seabass.gsfc.nasa.gov/wiki/FCHECK>
%    Just email fcheck@seabass.gsfc.nasa.gov with attached sb file(s).
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Basics
if nargin == 4
  num_data_types = 2; % full resolution (data_proc) and gridded (data_grid)
else
  num_data_types = 1; % Just full resolution (data_proc)
end

seabass_fields = {'PSD_DNSD' 'PSD_DVSD' 'VSF'};
ffmt = '%.4f'; % basic formatting string

%% 1 | Generate filename for writing
if isfield(cfg.header,'sb_filename')
  w.filename = cfg.header.sb_filename;
elseif isfield(cfg.header,'cruise_qualifer')
  w.filename = [cfg.header.experiment '-' cfg.header.cruise '_' cfg.header.instrument_model '_' cfg.header.cruise_qualifer '_' datestr(min(data_proc.date),'yyyymmdd') '-' datestr(max(data_proc.date),'yyyymmdd') '_R' cfg.header.release_number '.sb'];
else
  w.filename = [cfg.header.experiment '-' cfg.header.cruise '_' cfg.header.instrument_model '_' datestr(min(data_proc.date),'yyyymmdd') '-' datestr(max(data_proc.date),'yyyymmdd') '_R' cfg.header.release_number '.sb'];
end

% Pull out calibration files
calfiles = unique(data_proc.zscfile);

%% 2 | Loop through full resolution and gridded data to write sb formatted files
for dtype = 1:num_data_types
  
  %%
  w.fields_list = {'date'     'time'    'station' 'lat'     'lon'     'depth' 'c'   'trans' 'quality'};
  w.units_list  = {'yyyymmdd' 'HH:MM:SS' 'none'   'degrees' 'degrees' 'm'     '1/m' '%'     'none'};
  w.writefmt    = {'%s'       '%s'       '%d'     ffmt      ffmt      ffmt     ffmt '%.1f'  '%d'};
  
  %% 3 | Define lisst structure and populate basic fields
  lisst = struct();
  if dtype == 1 % FULL RESOLUTION
    data = data_proc;
  elseif dtype == 2 % GRIDDED
    %fprintf('Gridding may not be working...\n');
    data = data_grid;
  end
  
  % Cannot use data.date(:) because matlab assumes columnwise indexing,
  % must reshape the transpose since our data is organized in [cast,depth]
  
  date_array = reshape(data.date',[],1);
  try
    lisst.date = cellstr(datestr(date_array,'yyyymmdd'));
    lisst.time = cellstr(datestr(date_array,'HH:MM:SS'));
  catch
    fprintf('times may be nans -- indicating further problems...\n');
    keyboard
  end
  if dtype == 1 % FULL RESOLUTION
    lisst.station = data.cast;
    lisst.lat     = data.lat;
    lisst.lon     = data.lon;
    lisst.depth   = data.depth;
    lisst.c       = data.beam_attenuation; % Beam attenuation coefficient (NOT Beam attenuation coefficient of particles (ap + bp, e.g Absorption coefficient of particles (ad + aph) + Scattering coefficient of particles))
    lisst.trans   = data.tau.*100;         % Percent transmission [%]
    lisst.quality = data.quality_flag;     % Analyst-specific data quality flag. A definition of the quality flags should be provided as metadata header comments and within accompanying documentation files.
  elseif dtype == 2 % GRIDDED
    num_casts  = size(data.date,1);
    num_depths = size(data.date,2);
    %lisst.datenum = repelem(data.datenum, num_depths);
    lisst.station = repelem(data.cast, num_depths);
    lisst.lat     = repelem(data.lat,  num_depths);
    lisst.lon     = repelem(data.lon,  num_depths);
    lisst.depth   = repmat(data.depth',num_casts,1);
    lisst.c       = reshape(data.beam_attenuation',  [],1); % Beam attenuation coefficient (NOT Beam attenuation coefficient of particles (ap + bp, e.g Absorption coefficient of particles (ad + aph) + Scattering coefficient of particles))
    lisst.trans   = reshape((data.tau.*100)',[],1);         % Percent transmission [%]
    % Remove quality_flag because any failed points are ignored during
    % gridding
    %cfg.grid_options.ignore_flagged = 1 = throws out flagged data before gridding,
    %cfg.grid_options.ignore_flagged = 0 = grids all data regardless if flagged or not
    if cfg.grid_options.ignore_flagged
      %fprintf('STOP AND SEE WHAT QUALITY FLAGS LOOK LIKE!\n')
      %fprintf(' i.e. data was binned/gridded using nanmean, so flags may not be integers anymore...\n')
      lisst.quality = reshape(data.quality_flag',[],1);
      lisst.quality = round(lisst.quality);
      lisst.quality(isnan(lisst.quality)) = 2;
    else
      idx_remove = strcmp(w.fields_list,'quality');
      w.fields_list(idx_remove) = [];
      w.units_list(idx_remove)  = [];
      w.writefmt(idx_remove)    = [];
    end
    
    % add bincount where bincount = Number of records averaged into a bin or reported measurement
    lisst.bincount       = reshape(data.bincount',[],1);
    w.fields_list(end+1) = {'bincount'};
    w.units_list(end+1)  = {'none'};
    w.writefmt(end+1)    = {'%d'};
    % update sb filename to reflect binning
    w.filename = insertAfter(w.filename,'LISST-Deep','_binned');
    
    % Add to comments about bin depth
    cfg.header.comments = [cfg.header.comments; ['! depth represents the center of ' num2str(cfg.grid_options.bin_depth_m) ' meter depth bins']; '!'];
  end
  
  % add r2r_event, if possible
  if isfield(data,'r2r_event')
    if dtype == 1; lisst.R2R_Event = data.r2r_event; end
    if dtype == 2; lisst.R2R_Event = repelem(data.r2r_event, num_depths, 1); end
    w.fields_list = [w.fields_list 'R2R_Event'];
    w.units_list  = [w.units_list  'none'];
    w.writefmt    = [w.writefmt    '%s'];
  end
  
  %% 4 | Loop through SeaBASS fields and format for writing
  for sf = 1:numel(seabass_fields)
    switch seabass_fields{sf}
      case 'VSF'
        % VSF [1/m/sr]	Total volume scattering function (VSF### or
        % VSF###_YYYang, where YYY is collection angle). Also called beta.
        % generate fieldnames
        vsf_fields = cellstr(strcat('VSF',num2str(cfg.inst.wavelength),'_',num2str(meta_proc.VSF.angles.(['type' cfg.inst.typeX '_Median'])),'ang'))';
        vsf_fields = erase(vsf_fields,' ');
        % generate units
        if ~strcmp(meta_proc.VSF.unit,'1/m/sr')
          error('VSF needs units of 1/m/sr but unexpected units found so must do conversion\n')
        end
        vsf_units  = repmat({meta_proc.VSF.unit},size(vsf_fields));
        w.fields_list = [w.fields_list vsf_fields];
        w.units_list  = [w.units_list  vsf_units];
        w.writefmt    = [w.writefmt    repmat({ffmt},size(vsf_fields))];
        % Loop through size bins and allocate data to lisst structure
        for ns = 1:32
          if dtype == 1; tmp = data.VSF(:,ns);   end % FULL RESOLUTION
          if dtype == 2; tmp = data.VSF(:,:,ns); end % GRID RESOLUTION
          lisst.(strrep(vsf_fields{ns},'.','_')) = reshape(tmp',[],1);
        end
      case 'PSD_DNSD'
        % PSD_DNSD [number/m^3/um]	Differential number size distribution. Must
        % be used in conjunction with a suffix to specify sizes, e.g.,
        % PSD_DNSD_80.5umsize (see _###umsize in table of field name modifiers).
        % Dataset should include the /PSD_bin_size_method and
        % /PSD_bin_size_boundaries metadata headers.
        % generate fieldnames
        nsd_fields = cellstr(strcat('PSD_DNSD_',num2str(cfg.inst.bins.mid_point,'%.2f'),'umsize'))';
        nsd_fields = erase(nsd_fields,' ');
        % generate units
        if ~strcmp(meta_proc.PSD_DNSD.unit,'#/m^3/um') && ~strcmp(meta_proc.PSD_DNSD.unit,'number/m^3/um')
          error('PSD_DNSD needs units of number/m^3/um but unexpected units found so must do conversion\n')
        end
        w.fields_list = [w.fields_list nsd_fields];
        %w.units_list  = [w.units_list  repmat({meta_proc.PSD_DNSD.unit},size(nsd_fields))];
        w.units_list  = [w.units_list  repmat({'number/m^3/um'},size(nsd_fields))];
        w.writefmt    = [w.writefmt    repmat({ffmt},size(nsd_fields))];
        % Loop through size bins and allocate data to lisst structure
        for ns = 1:32
          if dtype == 1; tmp = data.PSD_DNSD(:,ns);   end % FULL RESOLUTION
          if dtype == 2; tmp = data.PSD_DNSD(:,:,ns); end % GRID RESOLUTION
          lisst.(strrep(nsd_fields{ns},'.','_')) = reshape(tmp',[],1);
        end
      case 'PSD_DVSD'
        % PSD_DVSD [ul/m^3/um]	Differential volume size distribution. Must be
        % used in conjunction with a suffix to specify sizes, e.g.,
        % PSD_DVSD_80.5umsize (see _###umsize in table of field name
        % modifiers).  Dataset should include the /PSD_bin_size_method and
        % /PSD_bin_size_boundaries metadata headers.
        % generate fieldnames
        vsd_fields = cellstr(strcat('PSD_DVSD_',num2str(cfg.inst.bins.mid_point,'%.2f'),'umsize'))';
        vsd_fields = erase(vsd_fields,' ');
        % generate units
        if ~strcmp(meta_proc.PSD_DVSD.unit,'ul/m^3/um') && ~strcmp(meta_proc.PSD_DVSD.unit,'uL/m^3/um')
          error('PSD_DVSD needs units of ul/m^3/um but unexpected units found so must do conversion\n')
        end
        w.fields_list = [w.fields_list vsd_fields];
        w.units_list  = [w.units_list  repmat({meta_proc.PSD_DVSD.unit},size(vsd_fields))];
        w.writefmt    = [w.writefmt    repmat({ffmt},size(vsd_fields))];
        % Loop through size bins and allocate data to lisst structure
        for ns = 1:32
          if dtype == 1; tmp = data.PSD_DVSD(:,ns);   end % FULL RESOLUTION
          if dtype == 2; tmp = data.PSD_DVSD(:,:,ns); end % GRID RESOLUTION
          lisst.(strrep(vsd_fields{ns},'.','_')) = reshape(tmp',[],1);
        end
      otherwise
        fprintf('Need to set this up! %s\n',seabass_fields{sf})
        keyboard
    end % HANDLE FIELD
  end % LOOP THROUGH SEABASS FIELDS TO WRITE
  
  %% 5 | Remove erroneous entires where date is not finite - i.e. no real data in that depth bin
  bad = isnan(date_array);
  fnm = fieldnames(lisst);
  for nf = 1:numel(fnm)
    lisst.(fnm{nf})(bad) = [];
  end
  
  %%  Add quality flag definitions to header comments
  if isfield(lisst,'quality')
    % Example from CTD sb file
    % ! Flags: 0 = good; 1 = flagged bad by wildedit or loopedit
    % Generate string to add to comments
    flags = unique(lisst.quality);
    
    for nflag = 1:numel(flags)
      iflag = meta_proc.quality_flag.flag_values == flags(nflag);
      if nflag == 1
        flag_str = ['! Flags: ' num2str(meta_proc.quality_flag.flag_values(iflag)) ' = ' meta_proc.quality_flag.flag_meanings{iflag}];
      else
        flag_str = [flag_str '; ' num2str(meta_proc.quality_flag.flag_values(iflag)) ' = ' meta_proc.quality_flag.flag_meanings{iflag}];
      end
    end
    if cfg.qaqc_options.expert_qc == 0
      flag_str = {'! Flags set using automatic QC tests recommended by Sequoia Scientific based on transmission and laser power. ';...
        '! Since only automatic QC tests were performed, no other quality control or quality assurance is given and quality is set to 2 (not evaluated)';...
        flag_str};
    end
    cfg.header.comments = [cfg.header.comments; flag_str; '!'];
  end
  
  %% 6 | Generate SeaBASS header text
  % pull out max & min date, time
  idx_good = isfinite(data.date(:));
  [~,imin] = nanmin(data.date(idx_good));
  [~,imax] = nanmax(data.date(idx_good));
  
  % Write formatted header for sb file
  hdr_SEABASS={'/begin_header';
    ['/investigators='      cfg.header.investigators];
    ['/affiliations='       cfg.header.affiliations];
    ['/contact='            cfg.header.contact];
    ['/experiment='         cfg.header.experiment];
    ['/cruise='             cfg.header.cruise];
    ['/data_file_name='     w.filename];
    ['/documents='          strjoin(cfg.header.documents,',')];
    ['/data_type='          cfg.header.data_type];
    ['/data_status='        cfg.header.data_status];
    ['/start_date='         lisst.date{imin}];
    ['/end_date='           lisst.date{imax}];
    ['/start_time='         lisst.time{imin} '[GMT]'];
    ['/end_time='           lisst.time{imax} '[GMT]'];
    ['/north_latitude='     num2str(max(data.lat),'%.4f') '[DEG]'];
    ['/south_latitude='     num2str(min(data.lat),'%.4f') '[DEG]'];
    ['/east_longitude='     num2str(max(data.lon),'%.4f') '[DEG]'];
    ['/west_longitude='     num2str(min(data.lon),'%.4f') '[DEG]'];
    '/water_depth=NA';
    '/station=NA';
    ['/missing='                  cfg.header.missing];
    ['/delimiter='                cfg.header.delimiter];
    ['/instrument_manufacturer='  cfg.header.instrument_manufacturer];
    ['/instrument_model='         cfg.header.instrument_model];
    ['/calibration_files='        strjoin(calfiles,',')];
    ['/calibration_date='         cfg.header.caldates];
    ['/PSD_bin_size_method='      strrep(cfg.header.PSD_bin_size_method,' ','_')]; %specify method used to select the nominal bin size such as arithmetic_mean or geometric_mean, just make sure there are no spaces in the string
    ['/PSD_bin_size_boundaries='  cfg.header.PSD_bin_size_boundaries];  %please provide a comma-separated list with the bin size boundaries in increasing order, e.g., 5,9.5,15,20 ...
    '! The sizes reported in the field names represent the bin center in micrometers and defined in PSD_bin_size_method.'
    '! size of the bin limits in micrometers';
    '!'};
  % Insert comments, then finish with /fields and /units
  hdr_SEABASS = [hdr_SEABASS; cfg.header.comments;
    
  '!';
  ['/fields=' strjoin(w.fields_list,',')];
  ['/units=' strjoin(w.units_list,',')];
  '/end_header'];

% check if there is whitespace in any metadata headers
% whitespace in comments is OKAY
if any(contains(hdr_SEABASS,' ') & ~contains(hdr_SEABASS,'!'))
  fprintf('White space in metadata header, must remove to pass fcheck\n')
  keyboard
end

%% 7 | Format for writing to sb file
w.writefmt = strjoin(w.writefmt,',');
w.writefmt = [w.writefmt '\n']; % add end of line character
% reformat odv2 table to temporary variable to enable simply write to file
dat_table = struct2table(lisst);
% Renove erroneous depths
dat_table(dat_table.depth < 0,:) = [];
if dtype == 1
  %dat_table = dat_table(1:5000,:);
elseif dtype == 2
  dat_table(dat_table.bincount == 0,:) = [];
end
% Get rid of erroneous data lines...
iPSD = contains(w.fields_list,'PSD');
PSDs = table2array(dat_table(:,iPSD));
ibad = all(isnan(PSDs),2);
dat_table(ibad,:) = [];

dat_write = table2cell(dat_table); % convert to cell array to handle different variable types
dat_write = dat_write';            % transpose because fprintf function prints data columnwise
% convert NaNs to missing identifier
dat_write(cellfun(@(x) isnumeric(x) && isnan(x), dat_write)) = {cfg.header.missing}; % missing=-9999 SeaBASS required header field

%% 8 | Write data to file
fprintf('\n  Writing data to: %s\n', fullfile(cfg.path.dir_submit_L2,w.filename));
fileID = fopen(fullfile(cfg.path.dir_submit_L2,w.filename),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile(cfg.path.dir_submit_L2,w.filename))
  keyboard
end
fprintf(fileID,'%s\n',hdr_SEABASS{:});  % write header
fprintf(fileID,w.writefmt,dat_write{:});     % write data
fclose(fileID);                 % close file

%   fprintf('\n  Writing data to screen\n');
%   fprintf('%s\n',hdr_SEABASS{:});  % write header
%   fprintf(w.writefmt,dat_write{:});     % write data

end % LOOP THROUGH DIFFERENT KINDS OF DATA - FULL RESOLUTION VS GRIDDED



end %% MAIN FUNCTION