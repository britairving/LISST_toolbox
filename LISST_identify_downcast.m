function [cfg, data_pre, meta_pre] = LISST_identify_downcast(cfg,data_pre,meta_pre,ctd)
%FUNCTION LISST_identify_downcast
%
%  Syntax:
%    [cfg, data_pre, meta_pre] = LISST_identify_downcast(cfg,data_pre,meta_pre,ctd)
%
%  Description:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
dbstop if error
ask_to_verify_downcasts_from_file = 1; % Default to yes
if strcmp(datestr(now,'yyyymmdd'),'20210408')
  ask_to_verify_downcasts_from_file = 0;
end
%% Set up
done_selecting_cfg.proc_options.profile_limit = 0;
while ~done_selecting_cfg.proc_options.profile_limit
  fprintf('Finding downcast section of each profile: current profile limit set to %d\n',cfg.proc_options.profile_limit)
  plim_chc = input(' Is that correct? <1/0> ');
  if plim_chc == 0
    cfg.proc_options.profile_limit = input(' Enter profile limit: ');
  else
    done_selecting_cfg.proc_options.profile_limit = 1;
  end
end

%% Initialize new variables
data_pre.profile_index = zeros(size(data_pre.date));
% Update description to reflect corrections
meta_pre.profile_index.name = 'logical indicies of downcast';
meta_pre.profile_index.values  = [1 0];
meta_pre.profile_index.meaning = {'1=downcast' '0=upcast or soak'};

%% 1 | Open profile index file for writing
% file that indicies for start and end of downcast for each cast
if exist(cfg.path.file_downcast,'file')
  downcasts = readtable(cfg.path.file_downcast);
  fileID = fopen(cfg.path.file_downcast,'a');
else
  downcasts.cast = nan(size(cfg.proc_options.cast));
  fileID = fopen(cfg.path.file_downcast,'a');
  fprintf(fileID,'cast,datfile,length_cast,downcast_start,downcast_end\n');
end

%% 1 | LOOP THROUGH LISST CASTS
% Skips this if all casts have been identified
makefig; ax1 = gca;
% Loop through casts
for nf = 1:numel(cfg.proc_options.cast)
  % reset switches
  find_downcast = 1;
  % Pull out cast number and associated lisst filename
  datname = cfg.proc_options.datfile{nf};
  ncast   = cfg.proc_options.cast(nf);
  scast   = num2str(ncast,'%03d');
  
  fprintf('\n-------------------------------------------------------------\n')
  fprintf('Finding downcast for cast#%d\n',ncast)
  % Pull out lisst data for this cast
  idx_lisst   = find(data_pre.cast == ncast);
  lisst_depth = data_pre.depth(idx_lisst);
  lisst_time  = data_pre.date(idx_lisst);
  % pull out ctd data for this cast
  idx_ctd   = find(ctd.cast == ncast);
  ctd_date  = ctd.date(idx_ctd);
  ctd_depth = ctd.depth(idx_ctd);
  
  %% filter depth
  filter_window = 45;
  smo_lisst = smooth(lisst_depth,filter_window);
  smo_ctd   = smooth(ctd_depth,filter_window);
  
  %% plot cast
  hlisst = plot(ax1,lisst_time,smo_lisst,'-','Color','k','LineWidth',3,'DisplayName','LISST profile');
  datetick(ax1,'x','HH:MM','keepticks'); grid(ax1,'on'); hold(ax1,'on'); ax1.YDir = 'reverse'; ax1.XTickLabelRotation = 45;
  title(ax1,{[datname ': ' datestr(data_pre.date(idx_lisst(1)),'mm/dd HH:MM') ' - ' datestr(data_pre.date(idx_lisst(end)),'mm/dd HH:MM') ];...
    ['Cast' scast ]});
  
  ylabel(ax1,'Depth [m]');
  min_x = min([lisst_time; ctd_date]);
  max_x = max([lisst_time; ctd_date]);
  xlim(ax1,[min_x-0.001 max_x+0.001]);
  try
    hlegend = legend(ax1,'show','Location','se'); hlegend.FontSize = 14;
  catch
    htext2 = text(ax1,0.9,0.10,'- LISST','Color','k','units','normalized','FontWeight','bold','BackgroundColor','w');
  end
  if strcmp(cfg.project,'LISST_sn4025_2017_ASGARD_SKQ201709S')
    xlim auto;
  end
  
  %% FIND DOWNCAST - EITHER BY LOADING, AUTOMATIC DETECTION, OR MANUALLY
  %% 1 | Load previous entries
  if ismember(ncast,downcasts.cast)
    ic = find(downcasts.cast == ncast); % find index of matching cast
    if ~isequal(downcasts.length_cast(ic),length(idx_lisst))
      fprintf('length of cast #%d does not match... remove from file: %s\n',ncast,cfg.path.file_downcast)
      find_downcast = 1;
    else
      if numel(ic) > 1 && all(isequal(downcasts.cast(ic),downcasts.cast(ic(1))))
        ic = ic(1);
      elseif numel(ic) > 1 && ~all(isequal(downcasts.cast(ic),downcasts.cast(ic(1))))
        fprintf('multiple matches in the downcast indice file for cast: %d... \n',ncast)
        keyboard
      end
      find_downcast = 0;
    end
    if find_downcast == 0 && ask_to_verify_downcasts_from_file == 0
      idx_downcast = idx_lisst(downcasts.downcast_start(ic):downcasts.downcast_end(ic));
      data_pre.profile_index(idx_downcast) = 1;
    elseif find_downcast == 0 
      try
        idx_downcast = idx_lisst(downcasts.downcast_start(ic):downcasts.downcast_end(ic));
        hp = plot(ax1,data_pre.date(idx_downcast),data_pre.depth(idx_downcast),'m-','LineWidth',3,'DisplayName','LISST downcast');
        pchc = input('Is the highlighted profile correct? <1/0> ');
        if pchc == 0 % INCORRECT
          find_downcast = 1;
          delete(hp);
        else % PROFILE IS CORRECT, DEFAULT
          find_downcast = 0;
          data_pre.profile_index(idx_downcast) = 1;
        end
      catch
        find_downcast = 1;
        delete(hp)
      end
    end
  end
  %% FIND DOWNCAST
  if find_downcast
    %% Try to identify profile automatically
    % First, find approximate downcast based on CTD. CTD files should not be
    % fully processed (i.e. binned and only downcast) but they likely have
    % already removed the soaking period, so base limit based on that
    % assumption.
    try
      idx_afterctd = find(lisst_time >= ctd_date(1),1);
      [~,imax] = max(smo_lisst);
      % Catch if stays at bottom for more than a few minutes - if so, try and
      % identify when ctd hits the deepest depth
      max_depth = smo_lisst(imax);
      idx_atdeepest = find(smo_lisst >= max_depth-2 & smo_lisst <= max_depth);
      min_atdeepest = etime(datevec(lisst_time(idx_atdeepest(end))), datevec(lisst_time(idx_atdeepest(1))))/60; % seconds to minutes
      if min_atdeepest > 1
        imax = find(smo_lisst > max_depth-1,1,'first');
      end
    catch
      keyboard
    end
    idx_approximate_downcast = idx_afterctd:imax;
    % Returns 0.5 for turning point, 1 for downcast, 0 for upcast
    try
      [profile_index, ~] = bi_find_downcast(smo_lisst(idx_approximate_downcast),lisst_time(idx_approximate_downcast), cfg.proc_options.profile_limit, 0);
      profile_index(profile_index ~= 1) = 0; % 1 = downcast, 0 = otherwise
    catch
      profile_index = zeros(size(idx_approximate_downcast));
    end

    if all(profile_index == 0)
      idx_downcast = idx_approximate_downcast;
    else
      idx_downcast = idx_approximate_downcast(profile_index == 1);
    end
    
    try
      hp = plot(ax1,lisst_time(idx_downcast),smo_lisst(idx_downcast),'m-','LineWidth',3,'DisplayName','LISST downcast');
    catch
      hp = plot(ax1,lisst_time,smo_lisst,'m-','LineWidth',3,'DisplayName','LISST downcast');
    end
    
    %% Decide if downcast is correct. If not, select manually
    pchc = input('Is the highlighted profile correct? <1/0> ');
    if isempty(pchc) || pchc == 1 % Default to YES if user hits <enter>
      data_pre.profile_index(idx_lisst(idx_downcast)) = 1;
    elseif pchc == 9
      keyboard
    else
      delete(hp);
      fprintf('\nSelect profile manually...\n')
      % Returns 1 for downcast, NaN otherwise
      [profile_index, ~] = select_downcast(lisst_time,smo_lisst);
      profile_index(profile_index ~= 1) = 0; % 1 = downcast, 0 = otherwise
      idx_downcast = profile_index == 1;
      hp = plot(ax1,lisst_time(idx_downcast),smo_lisst(idx_downcast),'m-','LineWidth',3,'DisplayName','LISST downcast');
      data_pre.profile_index(idx_lisst(idx_downcast)) = 1;
    end
    %% Write to file
    %fprintf(fileID,'cast,datfile,length_cast,downcast_start,downcast_end\n');
    if islogical(idx_downcast)
      downcast_start = find(idx_downcast,1,'first');
      downcast_end   = find(idx_downcast,1,'last');
    else
      downcast_start = idx_downcast(1);
      downcast_end   = idx_downcast(end);
    end
    fprintf(fileID,'%d,%s,%d,%d,%d\n',ncast,datname,length(idx_lisst),downcast_start,downcast_end);
    
  end %% FIND_DOWNCAST
  % save figure
  if cfg.savefig
    figname = fullfile(cfg.path.dir_figs,[cfg.project '_cast' scast '_downcast']);
    standard_printfig_lowrespng(figname)
  end
  cla; ax1 = gca;
end % LOOP THROUGH CASTS
close
end %% MAIN FUNCTION