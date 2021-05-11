function LISST_write_NGALTER(cfg,data_proc,meta_proc)
%FUNCTION LISST_write_NGALTER
%
%  Syntax:
%    LISST_write_NGALTER(cfg,data_grid,meta)
%  
%  Description:
%
%  Refereces:
%   https://docs.google.com/document/d/1ySx3GuoWYnDzqk6s2EVuytbooW_Y0u9DUYzEujtmiSg/edit#heading=h.wb34a31ccr3k
%
%  Notes:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Generate filename for writing
% I | NGA LTER Format | 'nga_lter' % AXIOM Research Workspace
wname = ['NGA_' cfg.header.cruise_ID '_LISST_L2_v1.csv'];
wname = fullfile(cfg.path.dir_submit_L2,wname);

%%
sz = size(data_proc.quality_flag);
ffmt = '%.4f'; % basic formatting string

%% 1 | Set up standard columns
hdr   = {'Cruise' 'Station' 'Type' 'Date_Time' 'Longitude_[decimal_degrees_east]' 'Latitude_[decimal_degrees_north]' 'Bot. Depth [m]' 'Cast_Number' 'Depth_[m]'};
units = {'none'   'none'    'none' 'yyyy-mm-ddTHH:MM:SS' 'degrees' 'degrees' 'm' 'none' 'm'};
wfmt  = {'%s'     '%s'      '%s'   '%s'                   ffmt      ffmt     '%d' '%d'  ffmt};
lisst = struct();
lisst.Cruise      = repmat({cfg.header.cruise_ID},sz);  %Sikuliaq standard or [vessel-code][YYYY][MM]
lisst.Station     = data_proc.station;                  %StationName
lisst.Type        = repmat({'C'},sz);                   %B for bottle, C for CTD, U for underway, etc
lisst.Date_Time   = cellstr(datestr(data_proc.date,'yyyy-mm-ddTHH:MM:SS'));  %Date and time of station in extended ISO8601 format ? 2008-04-23T15:15:00.000
lisst.Longitude   = data_proc.lon;                      %Longitude of station or measurement
lisst.Latitude    = data_proc.lat;                      %Latitude of station or measurement
lisst.Bot_depth   = data_proc.botdepth;                 %Bottom depth of station
lisst.Cast_Number = data_proc.cast;                   %Consecutive cast number
lisst.Depth       = data_proc.depth;                  %Depth of measurement, positive down, could be Pressure_[dbar]

%% 2 | Add LISST parameters
if isfield(meta_proc.PSD,'bins')
  bin_size = meta_proc.PSD.bins;
else
  fprintf('need to add bin_size here\n')
end
  
lisst_fields = {'tau' 'quality_flag' 'mean_size' 'PSD' 'VSD'};
for nl = 1:numel(lisst_fields)
  f = lisst_fields{nl};
  if strcmp(f,'tau')
    hdr = [hdr 'transmission_[%]'];  % Add name
    units = [units '%'];             % Add unit
    wfmt  = [wfmt '%.1f'];           % Add write format
    lisst.(f) = data_proc.(f)*100.0; % Add data
  elseif strcmp(f,'quality_flag')
    hdr = [hdr f];             % Add name
    units = [units 'none'];    % Add unit
    wfmt  = [wfmt '%d'];       % Add write format
    lisst.(f) = data_proc.(f); % Add data
  elseif size(data_proc.(f),2) == 1
    str = [f '_[' meta_proc.(f).unit ']']; 
    hdr = [hdr str];                    % Add name
    units = [units meta_proc.(f).unit]; % Add unit
    wfmt  = [wfmt ffmt];                % Add write format
    lisst.(f) = data_proc.(f);          % Add data
  elseif size(data_proc.(f),2) == numel(bin_size.slimits)
    for ns = 1:numel(bin_size.slimits)
      str = [f '_' bin_size.slimits{ns} '_[' meta_proc.(f).unit ']'];
      hdr = [hdr str];                              % Add name
      units = [units meta_proc.(f).unit];           % Add unit
      wfmt  = [wfmt ffmt];                          % Add write format
      lisst.([f num2str(ns)]) = data_proc.(f)(:,ns);% Add data
    end
  else
    fprintf('Unexpected size!\n')
    keyboard
  end
end
  
if strcmp(cfg.project,'LISST_sn4025_2019_NGA_TGX201909')
  t = struct2table(lisst);
  first_idx = find(strcmp(t.Properties.VariableNames,'mean_size'));
  datarray  = table2array(t(:,first_idx:end));
  idx_bad   = all(isnan(datarray),2);
  t(idx_bad,:) = [];
  lisst = table2struct(t,'ToScalar',true);
end

%% Format header for writing
cols = strjoin(hdr,',');
unit = strjoin(units,',');
fmt  = strjoin(wfmt,',');
fmt  = [fmt '\n'];

% reformat data table to temporary variable to enable simply write to file
lisst_write = struct2table(lisst);     % convert to table
lisst_write = table2cell(lisst_write); % convert to cell array to handle different variable types
lisst_write = lisst_write';            % transpose because fprintf function prints data columnwise
% % convert NaNs to missing identifier
% lisst_write(cellfun(@(x) isnumeric(x) && isnan(x), lisst_write)) = {-9999}; % missing=-9999 SeaBASS required header field

%% Write odv table to NONDIFFERENTIAL file
fprintf('\n  Writing data to: %s\n',wname);
fileID = fopen(wname,'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',wname)
  keyboard
end

fprintf(fileID,'%s\n',cols); % write header
fprintf(fileID,'%s\n',unit); % write units
fprintf(fileID,fmt,lisst_write{:});      % write data
fclose(fileID);                 % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
% fprintf('%s\n',cols) % write header
% fprintf('%s\n',unit); % write units
% fprintf(fmt,lisst_write{:})      % write data

%% Plot
LISST_plot_NGALTER(wname)

end %% MAIN FUNCTION 