function [cfg, data_pre, meta_pre] = LISST_match_ctd_cast(cfg, data_pre, meta_pre, ctd)
%FUNCTION LISST_match_ctd_cast
%
%  Syntax:
%    [cfg, data_pre] = LISST_match_ctd_cast(cfg,data_pre)
%
%  Description:
%    Project specific since ctd data formatting is different.
%
%  Refereces:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
dbstop if error
close all

remove_these_data = {}; % cell array to store lisst data filenames (data_pre.datfile) that should be removed
%% 0 | Initialize new fields
data_pre.cast = nan(size(data_pre.date));          % sequential cast number
data_pre.lat  = nan(size(data_pre.date));          % latitude in degrees north
data_pre.lon  = nan(size(data_pre.date));          % longitude in degrees east [0-180]
if isfield(ctd,'r2r_event')
  data_pre.r2r_event      = repmat({''},size(data_pre.date)); % r2r_event unique id
  meta_pre.r2r_event.name = 'Rolling Deck to Repository (R2R) Unique sample event number';
  meta_pre.r2r_event.unit = 'none';
end
cfg.proc_options.cast = nan(size(cfg.proc_options.datfile));

% Add to metadata
meta_pre.cast.name = 'sequential cast number';
meta_pre.cast.unit = 'none';
meta_pre.lat.name = 'Sample latitude (decimal fractions; -90 to 90)';
meta_pre.lat.unit = 'degrees';
meta_pre.lon.name = 'Sample longitude (decimal fractions; -180 to 180)';
meta_pre.lon.unit = 'degrees';

if isfield(ctd,'station')
  data_pre.station = repmat({''},size(data_pre.date)); 
  meta_pre.station.name = 'station name';
  meta_pre.station.unit = 'none';
end
if isfield(ctd,'botdepth')
  data_pre.botdepth = nan(size(data_pre.date));  
  meta_pre.botdepth.name = 'Bottom depth of station';
  meta_pre.botdepth.unit = 'm';
end
  
%% 1 | Open ctd match file for writing
if exist(cfg.path.file_ctdmatch,'file')
  castmatch = readtable(cfg.path.file_ctdmatch);
  fileID = fopen(cfg.path.file_ctdmatch,'a');
else
  castmatch = table();
  castmatch.cast    = nan(numel(cfg.proc_options.datfile),1); % Sequential cast number
  castmatch.datfile = cfg.proc_options.datfile';              % filename of raw LISST datafile (.dat extension)
  fileID = fopen(cfg.path.file_ctdmatch,'a');
  fprintf(fileID,'cast,datfile\n');
end

%% 1 | LOOP THROUGH UNIQUE LISST PFOFILES
% Skips this if all casts have been identified
skip_files = {};
makefig; ax1 = subplot(1,2,1); ax2 = subplot(1,2,2);
% Loop through lisst data and find matching cast
for nf = 1:numel(cfg.proc_options.datfile)
  datname   = cfg.proc_options.datfile{nf};
  fprintf('\n-------------------------------------------------------------\n')
  fprintf('Looking at data from %s\n',datname)
  % Pull out indicies for this cast/data file
  idx_lisst = find(strcmp(data_pre.datfile, datname));
  
  %% find nearest ctd data
  % Don't repeat work
  if ismember(datname,castmatch.datfile)
    ic = find(ismember(castmatch.datfile,datname));
    if numel(ic) > 1 && all(isequal(castmatch.cast(ic),castmatch.cast(ic(1))))
      ic = ic(1);
    elseif numel(ic) > 1 && ~all(isequal(castmatch.cast(ic),castmatch.cast(ic(1))))
      fprintf('multiple casts match with this lisst data: %s... \n',datname)
      keyboard
    end
    rng_ctd = find(ctd.cast == castmatch.cast(ic));
    if isempty(rng_ctd)
      cast_selected = 0;
    else
      data_pre.cast(idx_lisst) = ctd.cast(rng_ctd(1));
      data_pre.lat(idx_lisst)  = ctd.lat(rng_ctd(1));
      data_pre.lon(idx_lisst)  = ctd.lon(rng_ctd(1));
      if isfield(data_pre,'r2r_event')
        data_pre.r2r_event(idx_lisst) = ctd.r2r_event(rng_ctd(1));
      end
      if isfield(data_pre,'station')
        data_pre.station(idx_lisst) = ctd.station(rng_ctd(1));
      end
      if isfield(data_pre,'botdepth')
        data_pre.botdepth(idx_lisst) = ctd.botdepth(rng_ctd(1));
      end
      cfg.proc_options.cast(nf) = ctd.cast(rng_ctd(1));
      cast_selected = 1;
    end
  else
    cast_selected = 0;
  end
  
  %% 1 | Find matching cast
  skip_matchctd  = 0;
  while ~cast_selected
    
    % Pull out time
    st1 = datestr(data_pre.date(idx_lisst(1)),'mm/dd HH:MM');
    st2 = datestr(data_pre.date(idx_lisst(end)),    'HH:MM');
    
    %% plot LISST data
    plot(ax1,data_pre.date(idx_lisst),data_pre.depth(idx_lisst),'b.','DisplayName',datname)
    datetick(ax1,'x','HH:MM','keepticks'); grid(ax1,'on'); hold(ax1,'on'); ax1.YDir = 'reverse';
    text(ax1,0.98,0.12,{'LISST profile';[st1 '-' st2]},'Color','b','units','normalized','HorizontalAlignment','right');
    ax1.XTickLabelRotation = 45;
    title(ax1,[datname ': ' datestr(data_pre.date(1),'mm/dd HH:MM') ' - ' datestr(data_pre.date(end),'mm/dd HH:MM') ]);
    axes(ax1);pause(1);
    ylabel(ax1,'Depth [m]');
    
    plot(ax2,data_pre.temp(idx_lisst),data_pre.depth(idx_lisst),'b.-','DisplayName',datname)
    hold(ax2,'on'); grid(ax2,'on'); ax2.YDir = 'reverse';
    xlabel(ax2,'temperature'); %
    %hl2 = legend(ax2,'show','Location','se'); hl2.FontSize = 12;
    ax2.YLim = ax1.YLim;
    title(ax2,'Temperature not expected to match exactly')
    if skip_matchctd
      idx_match = find(ctd.cast == ncast);
    else
      idx_match = nearest(ctd.date,data_pre.date(idx_lisst(1)));
    end
    ncast = ctd.cast(idx_match(1));
    nlat  = ctd.lat(idx_match(1));
    nlon  = ctd.lon(idx_match(1));
    
    % range of ctd data
    rng_ctd = find(ctd.cast == ncast);
    if isempty(rng_ctd)
      fprintf('\nCast %d not found!\n',ncast)
      ncast = ncast + 1;
      fprintf('Trying to match cast %d not found!\n',ncast)
      rng_ctd = find(ctd.cast == ncast);
    end
    %plot ctd data
    hctd = plot(ax1,ctd.date(rng_ctd),ctd.depth(rng_ctd),'k.-','LineWidth',2,'DisplayName', ['CTD cast' num2str(ctd.cast(rng_ctd(1)))]);
    ht2 = text(ax1,0.98,0.04,{['CTD cast ' num2str(ncast)];[datestr(ctd.date(rng_ctd(1)),'mm/dd HH:MM') ' - ' datestr(ctd.date(rng_ctd(end)),'HH:MM') ]},'units','normalized','HorizontalAlignment','right','Color','k');
    axis(ax1,'tight')
    
    % plot temperatures
    plot(ax2,ctd.temp(rng_ctd),ctd.depth(rng_ctd),'k.-','LineWidth',2,'DisplayName',['cast' num2str(ncast)])
    axis(ax2,'tight');
    ax2.YLim = ax1.YLim;
    %% now decide what to do
    cast_selected  = 1; % assume correct cast
    done_selecting = 0;
    skip_matchctd  = 0;
    while ~done_selecting
      fprintf('Select %s...\n',datname)
      fprintf('  <0>  Ignore LISST profile (e.g. soak)\n')
      fprintf('  <1>  cast %d correct\n',ncast)
      fprintf('  <2>  cast %d incorrect, try next cast\n',ncast)
      fprintf('  <3>  cast %d incorrect, try last cast\n',ncast)
      fprintf('  <4>  skip for now\n')
      fprintf('  <5>  enter cast number manually\n')
      fprintf('  <99> stop\n')
      chc = input('  Enter choice: ');
      if isempty(chc); chc = 1; end
      switch chc
        case 0 % Ignore LISST profile
          remove_these_data = [remove_these_data; datname];
          done_selecting    = 1;
          cast_selected     = 1;
        case 1 % Cast correct
          done_selecting = 1;
          cast_selected  = 1;
          data_pre.cast(idx_lisst) = ncast;
          data_pre.lat(idx_lisst)  = nlat;
          data_pre.lon(idx_lisst)  = nlon;
          if isfield(data_pre,'r2r_event')
            data_pre.r2r_event(:) = ctd.r2r_event(idx_match(1));
          end
          if isfield(data_pre,'station')
            data_pre.station(:) = ctd.station(idx_match(1));
          end
          if isfield(data_pre,'botdepth')
            data_pre.botdepth(:) = ctd.botdepth(idx_match(1));
          end
          cfg.proc_options.cast(nf) = ncast;
        case 2 % try last cast
          ncast = ncast + 1;
          skip_matchctd  = 1;
          done_selecting = 1;
          cast_selected  = 0;
          delete(hctd);delete(ht2);cla(ax2);
          title(ax1,'');
        case 3 % Try next cast
          ncast = ncast - 1;
          done_selecting = 1;
          skip_matchctd  = 1;
          cast_selected  = 0;
          delete(hctd);delete(ht2);cla(ax2);
          title(ax1,'');
        case 4 % skip for now
          done_selecting = 1;
          cast_selected  = 1;
          % flag to come back to later, unclear...
          skip_files = [skip_files; datname];
        case 5 % Manual
          fprintf('\n  Chose to enter cast number manually...\n')
          ncast = input('  Enter cast number: ');
          done_selecting = 1;
          idx_match = find(ctd.cast == ncast);
          cfg.proc_options.cast(nf) = ncast;
          data_pre.cast(idx_lisst) = ncast;
          data_pre.lat(idx_lisst) = ctd.lat(idx_match(1));
          data_pre.lon(idx_lisst) = ctd.lon(idx_match(1));
          if isfield(data_pre,'r2r_event')
            data_pre.r2r_event(:) = ctd.r2r_event(idx_match(1));
          end
          if isfield(data_pre,'station')
            data_pre.station(:) = ctd.station(idx_match(1));
          end
          if isfield(data_pre,'botdepth')
            data_pre.botdepth(:) = ctd.botdepth(idx_match(1));
          end
        case 99
          fprintf('enter "dbcont to continue\n')
          keyboard
      end % Switch choice
    end % while done_selecting
    
    % Write to file
    if ~ismember(datname, remove_these_data) && cast_selected == 1 && chc ~=4
      fprintf(fileID,'%d,%s\n',ncast,datname);
    end
    
  end % while cast_selected

  % clear axes
  if exist('chc','var')
    if cfg.savefig && ~isequal(chc,3)
      figname = fullfile(cfg.path.dir_figs,[cfg.project '_' erase(datname,'.DAT') '_castmatch']);
      standard_printfig_lowrespng(figname)
    end
    cla(ax1);cla(ax2);
    title(ax1,'');title(ax2,'');
    clearvars chc
  end
  
end
if ~isempty(skip_files)
  fprintf('some files were skipped... decide what to do here\n')
  fprintf('%s\n',skip_files{:})
  keyboard
end
%% Move bad files
if ~isempty(remove_these_data)
  fprintf('removed these!!\n')
  ignore_dir = fullfile(cfg.path.dir_raw,'_ignore');
  for nrm = 1:numel(remove_these_data)
    if exist(fullfile(cfg.path.dir_raw,remove_these_data{nrm}),'file')
      movefile(fullfile(cfg.path.dir_raw,remove_these_data{nrm}),fullfile(ignore_dir,remove_these_data{nrm}));
    elseif exist(fullfile(ignore_dir,remove_these_data{nrm}),'file')
      fprintf('%s already in "_ignore" folder\n',remove_these_data{nrm})
    else
      fprintf('could not locate %s...\n',remove_these_data{nrm})
      keyboard
    end
  end
end
fclose(fileID);
close

%% Get rid of unwatned data - usually files generated during soak period
if ~isempty(remove_these_data)
  idx_remove_data = ismember(data_pre.datfile,remove_these_data);
  fields = fieldnames(data_pre);
  for nf = 1:numel(fields)
    data_pre.(fields{nf})(idx_remove_data,:) = [];
  end
  idx_remove_these = ismember(cfg.proc_options.datfile,remove_these_data);
  cfg.proc_options.datfile(idx_remove_these) = [];
  cfg.proc_options.zscfile(idx_remove_these) = [];
  cfg.proc_options.cast(idx_remove_these)    = [];
end

end %% MAIN FUNCTION