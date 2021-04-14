function LISST_write_SeaBASS(cfg,data_proc,data_grid,meta_proc)
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
if nargin == 3
  num_data_types = 2; % full resolution (data_proc) and gridded (data_grid)
else
  num_data_types = 1; % Just full resolution (data_proc)
end

seabass_fields = {'PSD_DNSD' 'PSD_DVSD' 'VSF'};
ffmt = '%.4f'; % basic formatting string

w.fields_list = {'date'     'time'    'station' 'lat'     'lon'     'depth' 'c'};
w.units_list  = {'yyyymmdd' 'HH:MM:SS' 'none'   'degrees' 'degrees' 'm'     '1/m'};
w.writefmt    = {'%s'       '%s'       '%d'     ffmt      ffmt      ffmt     ffmt };

  
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
  %% 3 | Define lisst structure and populate basic fields
  lisst = struct();
  if dtype == 1 % FULL RESOLUTION
    data = data_proc;
  elseif dtype == 2 % GRIDDED
    data = data_grid;
  end
  
  lisst.date = cellstr(datestr(data.date(:),'yyyymmdd'));
  lisst.time = cellstr(datestr(data.date(:),'HH:MM:SS'));
  if dtype == 1 % FULL RESOLUTION
    lisst.station = data.cast;
    lisst.lat     = data.lat;
    lisst.lon     = data.lon;
    lisst.depth   = data.depth; 
  elseif dtype == 2 % GRIDDED
    lisst.station = repmat(data.cast(:), size(data.date,2),1);
    lisst.lat     = repmat(data.lat(:),  size(data.date,2),1);
    lisst.lon     = repmat(data.lon(:),  size(data.date,2),1);
    lisst.depth   = repmat(data.depth(:),size(data.date,1),1);
    % add bincount
    lisst.bincount  = data.bincount(:);         % Number of records averaged into a bin or reported measurement
    w.fields_list(end:end+1) = {'bincount' 'c'};
    w.units_list(end:end+1)  = {'none' '1/m'};
    w.writefmt(end:end+1)    = {'%d' ffmt};
  end
  lisst.c = data.beam_attenuation(:); % Beam attenuation coefficient (NOT Beam attenuation coefficient of particles (ap + bp, e.g Absorption coefficient of particles (ad + aph) + Scattering coefficient of particles))
  
  % add r2r_event, if possible
  if isfield(data,'r2r_event')
    if dtype == 1; lisst.R2R_Event = data.r2r_event(:); end
    if dtype == 2; lisst.R2R_Event = repmat(data.r2r_event(:), size(data.date,2),1); end
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
          lisst.(strrep(vsf_fields{ns},'.','_')) = tmp(:);
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
        w.units_list  = [w.units_list  repmat({meta_proc.PSD_DNSD.unit},size(nsd_fields))];
        w.writefmt    = [w.writefmt    repmat({ffmt},size(nsd_fields))];
        % Loop through size bins and allocate data to lisst structure
        for ns = 1:32
          if dtype == 1; tmp = data.PSD_DNSD(:,ns);   end % FULL RESOLUTION
          if dtype == 2; tmp = data.PSD_DNSD(:,:,ns); end % GRID RESOLUTION
          lisst.(strrep(nsd_fields{ns},'.','_')) = tmp(:);
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
          lisst.(strrep(vsd_fields{ns},'.','_')) = tmp(:);
        end
      otherwise
        fprintf('Need to set this up! %s\n',seabass_fields{sf})
        keyboard
    end % HANDLE FIELD
  end % LOOP THROUGH SEABASS FIELDS TO WRITE 

  
  %% 5 | Remove erroneous entires where date is not finite
  bad = isnan(data.date(:));
  fnm = fieldnames(lisst);
  for nf = 1:numel(fnm)
    lisst.(fnm{nf})(bad) = [];
  end
  
  %% 6 | Generate SeaBASS header text
  % pull out max & min date, time
  [~,imin] = nanmin(data.date(:));
  [~,imax] = nanmax(data.date(:));
  
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
  dat_write = table2cell(dat_table); % convert to cell array to handle different variable types
  dat_write = dat_write';            % transpose because fprintf function prints data columnwise
  % convert NaNs to missing identifier
  dat_write(cellfun(@(x) isnumeric(x) && isnan(x), dat_write)) = {cfg.header.missing}; % missing=-9999 SeaBASS required header field

  
  %% 8 | Write data to file
  fprintf('\n  Writing data to: %s\n', fullfile(cfg.path.submit_dir,w.filename));
  fileID = fopen(fullfile(cfg.path.submit_dir,w.filename),'w');  % open file
  if fileID < 0
    fprintf(' *** Error opening file %s\n',fullfile(cfg.path.submit_dir,w.filename))
    keyboard
  end
  fprintf(fileID,'%s\n',hdr_SEABASS{:});  % write header
  fprintf(fileID,fmt,dat_write{:});     % write data
  fclose(fileID);                 % close file
  keyboard
end % LOOP THROUGH DIFFERENT KINDS OF DATA - FULL RESOLUTION VS GRIDDED



end %% MAIN FUNCTION 