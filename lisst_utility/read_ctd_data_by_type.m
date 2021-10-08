function [ctd,ctd_filename] = read_ctd_data_by_type(ctd_type,ctd_directory)
%FUNCTION READ_CTD_DATA_BY_TYPE
%
%  Description:
%     Read in CTD data given file type and directory
%
%  Syntax:
%     ctd = read_ctd_data_by_type(ctd_type,ctd_directory)
%
%  Inputs:
%     ctd_type      | string containing type of ctd file
%     ctd_directory | path to directory where ctd files are
%
%  Outputs:
%     ctd = structure with the following fields:
%     "station"   | Station Name
%     "cast"      | Unique identifier for cast, usually sequential
%     "lat"       | Latitude of CTD cast in decimal degrees
%     "lon"       | Longitude of CTD cast in decimal degrees
%     "date"      | Date of CTD cast in datenum format
%     "press"     | seawater pressure in dbar
%     "temp"      | Seawater temperature in degrees Celsius
%     "salt"      | Seawater salinity in PSU
%     "depth"     | Depth in meters (positive down)
%     "botdepth"  | Bottom depth of cast in meters (optional)
%     "r2r_event" | Rolling Deck to Repository (R2R) Unique sample event number
%     "file"      | Filename of matlab CTD data
%           
%  Authors:
%     Brita K Irving  <bkirving@alaska.edu>
%% 1 | Initialize filename
ctd_filename = fullfile(ctd_directory,['CTD_' ctd_type '.mat']);

keyboard
%% 2 | Read CTD data based on type
switch ctd_type
  case 'seabass'
    ctdfiles = dir(fullfile(ctd_directory,'*.sb'));
    bad = contains({ctdfiles.name},'._');
    ctdfiles(bad) = [];
    for nsb = 1:numel(ctdfiles)
      [data, sbHeader, ~]  = readsb(fullfile(ctd_directory,ctdfiles(nsb).name),'MakeStructure',1);
      if nsb == 1 % Initialize ctd structure
        ctd.cast    = repmat(sbHeader.station,size(data.date));
        ctd.lat     = repmat(sbHeader.north_latitude,size(data.date));
        ctd.lon     = repmat(sbHeader.east_longitude,size(data.date));
        ctd.file    = repmat({ctdfiles(nsb).name},size(data.date));
        if isfield(sbHeader,'r2r_event') || isfield(sbHeader,'R2R_Event') 
          try ctd.r2r_event = repmat({sbHeader.r2r_event},size(data.date));
          catch; ctd.r2r_event = repmat({sbHeader.R2R_Event},size(data.date));end
        end
        if any(contains(sbHeader.comments,'Station ID:'))
          station_str = sbHeader.comments(contains(sbHeader.comments,'Station ID'));
          station_str = strsplit(char(station_str),'Station ID:');
          ctd.station = repmat({strtrim(station_str{2})},size(data.date));
        end
        ctd.date   = data.datenum;  % datenum calculated in readsb.m
        ctd.press  = data.pressure; % water pressure [dbar]
        ctd.temp   = data.wt;       % water temperature [ degreesC]
        ctd.salt   = data.sal;      % salinity [PSU];
      else % Concatenate data to ctd structure
        ctd.cast    = [ctd.cast; repmat(sbHeader.station,size(data.date))];
        ctd.lat     = [ctd.lat;  repmat(sbHeader.north_latitude,size(data.date))];
        ctd.lon     = [ctd.lon;  repmat(sbHeader.east_longitude,size(data.date))];
        ctd.file    = [ctd.file; repmat({ctdfiles(nsb).name},size(data.date))];
        if isfield(ctd,'r2r_event') 
          try    ctd.r2r_event = [ctd.r2r_event; repmat({sbHeader.r2r_event},size(data.date))];
          catch; ctd.r2r_event = [ctd.r2r_event; repmat({sbHeader.R2R_Event},size(data.date))];end
        end
        if any(contains(sbHeader.comments,'Station ID:'))
          station_str = sbHeader.comments(contains(sbHeader.comments,'Station ID'));
          station_str = strsplit(char(station_str),'Station ID:');
          ctd.station = [ctd.station; repmat({strtrim(station_str{2})},size(data.date))];
        else
          ctd.station = [ctd.station; repmat({'unknown'},size(data.date))];
        end
        ctd.date   = [ctd.date;   data.datenum];  
        ctd.press  = [ctd.press;  data.pressure]; 
        ctd.temp   = [ctd.temp;   data.wt];       
        ctd.salt   = [ctd.salt;   data.sal];      
      end
     
    end
    %% Remove any duplicate entries
    [~,uidx] = unique(ctd.date); % This will also sort by date
    fields = fieldnames(ctd);
    for nf = 1:numel(fields)
      ctd.(fields{nf}) = ctd.(fields{nf})(uidx,:);
    end
    
  case 'nga_lter'
    fprintf('You will need to get raw hex files and process through SBE software\n')
    fprintf('It is important to have a continuous timestamp (timeJ) and upcast and downcast data\n')
    fprintf('Then, point to that directory in the NGALTER_read_ctd_L0.m script where cnv files are located\n')
    keyboard
  otherwise 
    fprintf('NEED TO SET THIS UP!\n')
    keyboard
end

%% 3 | Calculate depth from pressure and latitude
if ~isfield(ctd,'depth')
  ctd.depth = -1*gsw_z_from_p(ctd.press,ctd.lat);
end

%% 5 | Save CTD data to matlab file
fprintf('Saving CTD data to %s\n',ctd_filename)
save(ctd_filename,'ctd','-v7.3');
end %% MAIN FUNCTION
