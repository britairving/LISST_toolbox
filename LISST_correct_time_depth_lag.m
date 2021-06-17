function [cfg, data_pre, meta_pre] = LISST_correct_time_depth_lag(cfg,data_pre,meta_pre,ctd)
% %FUNCTION LISST_correct_time_depth_lag
% %
% %  Syntax:
% %    [cfg, data_pre] = LISST_correct_time_depth_lag(cfg,data_pre)
% %
% %  Description:
% %    Project specific since ctd data formatting is different.
% %
% %  Refereces:
% %  Jean-Pascal Rueff (2021). Click-n'-drag plot %
% %  (https://www.mathworks.com/matlabcentral/fileexchange/5751-click-n-drag-plot),
% %  MATLAB Central File Exchange. Retrieved June 14, 2021.
% %
% %  Authors:
% %    Brita K Irving  <bkirving@alaska.edu>
%
%% 0 | Set up base structure
close all
dbstop if error
cfg.savefig = 0;
auto_methods = 2;
%% Initialize new variables
data_pre.depth_orig = data_pre.depth;
data_pre.date_orig  = data_pre.date;
meta_pre.depth_orig = meta_pre.depth;
meta_pre.date_orig  = meta_pre.date;

% Update description to reflect corrections
meta_pre.depth.name      = 'depth corrected for linear drift';
meta_pre.date.name       = 'MATLAB datenum corrected for time lag';

%% 3 | Find time and/or depth lag
% do this in separate loop so don't hold up the cast identification process

if ~isfield(cfg.proc_options,'lag_correction')
  cfg.proc_options.cast       = cfg.proc_options.cast';           % sequential cast number
  cfg.proc_options.datfile    = cfg.proc_options.datfile';        % LISST profile ( datafile name )
  cfg.proc_options.time_lag   = nan(size(cfg.proc_options.cast)); % [seconds] different between CTD time (assumed to be correct) and LISST time in seconds
  cfg.proc_options.depth1     = nan(size(cfg.proc_options.cast)); % [meter]   LISST depth at the first depth
  cfg.proc_options.depth2     = nan(size(cfg.proc_options.cast)); % [meter]   LISST depth at the second depth
  cfg.proc_options.depth1_lag = nan(size(cfg.proc_options.cast)); % [meter]   different between CTD depth (assumed to be correct) and LISST depth at the first depth
  cfg.proc_options.depth2_lag = nan(size(cfg.proc_options.cast)); % [meter]   different between CTD depth (assumed to be correct) and LISST depth at the second depth
end

%% 1 | Open ctd match file for writing
% Hopefully will not have to repeat work because of this
if exist(cfg.path.file_ctdlag,'file')
  ctdlag = readtable(cfg.path.file_ctdlag);
  fileID = fopen(cfg.path.file_ctdlag,'a');
else
  ctdlag.cast = nan(size(cfg.proc_options.cast));
  fileID = fopen(cfg.path.file_ctdlag,'a');
  fprintf(fileID,'cast,datfile,time_lag,depth1,depth2,depth1_lag,depth2_lag\n');
end

% if strcmp(datestr(now,'yyyymmdd'),'20210407')
%   ctdlag.cast = [];
% end

%% 1 | LOOP THROUGH LISST CASTS
% Skips this if all casts have been identified
makefig; ax1 = gca;
% Loop through casts
for nf = 1:numel(cfg.proc_options.cast)

  % Reset switches
  skip_auto = 0;
  ncnt_auto = 1;
  select_manually = 0;
  
  % Pull out cast number and associated lisst filename
  datname = cfg.proc_options.datfile{nf};
  ncast   = cfg.proc_options.cast(nf);
  scast   = num2str(ncast,'%03d');
  
  fprintf('\n-------------------------------------------------------------\n')
  fprintf('Correction cast#%d data\n',ncast)
  % Pull out lisst data for this cast
  idx_lisst = find(data_pre.cast == ncast);
  lisst_depth = data_pre.depth(idx_lisst);
  lisst_time  = data_pre.date(idx_lisst);
  % pull out ctd data for this cast
  idx_ctd   = find(ctd.cast == ncast);
  ctd_date  = ctd.date(idx_ctd);
  ctd_depth = ctd.depth(idx_ctd);
  
  %% Interpolate nans in the lisst data
  % interpolate over any nans at the ends of the time vector if they exist
  try
    lisst_time = interp1( find( ~isnan(lisst_time)),  lisst_time(~isnan(lisst_time)), 1:length(lisst_time), 'linear', 'extrap')';
    % do the same thing with depth
    lisst_depth = interp1( find( ~isnan(lisst_depth)), lisst_depth(~isnan(lisst_depth)), 1:length(lisst_depth), 'linear', 'extrap')';
  catch
    keyboard
  end
  %% filter depth
  smo_lisst = smooth(lisst_depth,25); % smooth lisst data less because resolution will be less than CTD usually
  smo_ctd   = smooth(ctd_depth,25);
  
  %% plot LISST data
  hlisst = plot(ax1,lisst_time,smo_lisst,'-','Color','k','LineWidth',3,'DisplayName','LISST profile');
  datetick(ax1,'x','HH:MM','keepticks'); grid(ax1,'on'); hold(ax1,'on'); ax1.YDir = 'reverse';
  ax1.XTickLabelRotation = 45;
  % highlight lisst profile
  title(ax1,{[datname ': ' datestr(data_pre.date(idx_lisst(1)),'mm/dd HH:MM') ' - ' datestr(data_pre.date(idx_lisst(end)),'mm/dd HH:MM') ];...
    ['CTD Cast' scast ]});
  ht1 = text(ax1,0.98,0.12,{'LISST profile';[datestr(lisst_time(1),'mm/dd HH:MM') ' - ' datestr(lisst_time(end),'HH:MM') ]},'units','normalized','HorizontalAlignment','right','Color','k','FontWeight','bold');
  try
    hlegend = legend(ax1,'show','Location','se'); hlegend.FontSize = 14;
  end
  ylabel(ax1,'Depth [m]');
  min_x = min([lisst_time; ctd_date]);
  max_x = max([lisst_time; ctd_date]);
  xlim(ax1,[min_x-0.001 max_x+0.001]);
  % plot ctd
  hctd = plot(ax1,ctd_date,smo_ctd,'b.-','LineWidth',6,'DisplayName', ['CTD cast' num2str(scast)]);
  ht2 = text(ax1,0.98,0.04,{['CTD cast ' scast];[datestr(ctd_date(1),'mm/dd HH:MM') ' - ' datestr(ctd_date(end),'HH:MM') ]},'units','normalized','HorizontalAlignment','right','Color','b','FontWeight','bold');
  
  %% CORRECION ALREADY FOUND AND STORED
  if ismember(ncast,ctdlag.cast)
    ic = find(ctdlag.cast == ncast); % find index of matching cast
    if numel(ic) > 1 && all(isequal(ctdlag.cast(ic),ctdlag.cast(ic(1))))
      ic = ic(1);
    elseif numel(ic) > 1 && ~all(isequal(ctdlag.cast(ic),ctdlag.cast(ic(1))))
      fprintf('multiple casts match with this lisst data: %s... \n',datname)
      keyboard
    end
    % store values to structure
    time_lag_seconds = ctdlag.time_lag(ic);
    dlag_depth1      = ctdlag.depth1(ic);
    dlag_depth2      = ctdlag.depth2(ic);
    dlag1            = ctdlag.depth1_lag(ic);
    dlag2            = ctdlag.depth2_lag(ic);
    
    % Correct depth for lag assuming linear drift
    try
      depth_bin_estimate = 0:1:fix(max(smo_lisst));
      if depth_bin_estimate == 0
        depth_bin_estimate = 0:1;
      end
      depth_lag_estimate = interp1([dlag_depth1 dlag_depth2],[dlag1 dlag2],depth_bin_estimate,'linear','extrap');
      depth_lag = interp1(depth_bin_estimate,depth_lag_estimate,smo_lisst);
      depth_lag = fillmissing(depth_lag,'linear','EndValues','extrap');
      depth_cor =  smo_lisst + depth_lag;
    catch
      keyboard
    end
    % time lag - convert from seconds to datenum (days)
    time_difference = time_lag_seconds./60./60./24; % convert from seconds to datenum
    lisst_time_cor  = lisst_time + time_difference;
    
    h1_cor = plot(ax1,lisst_time_cor,depth_cor,'g.-','LineWidth',1,'DisplayName','LISST corrected');
    hlagt = text(ax1,0.02,0.15,['time offset  = ' num2str(time_lag_seconds,'%.2f') ' seconds'],'units','normalized','FontWeight','bold');
    hlag1 = text(ax1,0.02,0.10, ['depth offset = ' num2str(dlag1,'%.2f')  'm @ ' num2str(dlag_depth1,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
    hlag2 = text(ax1,0.02,0.05, ['depth offset = ' num2str(dlag2,'%.2f') 'm @ ' num2str(dlag_depth2,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
    % save figure
    if cfg.savefig
      figname = fullfile(cfg.path.dir_figs,[cfg.project '_cast' scast '_sensor_lag']);
      standard_printfig_lowrespng(figname)
    end
    
  else
    %% FIND CORRECTION
    if skip_auto
      auto_done = 1;
      select_manually = 1;
    else
      auto_done = 0;
    end
    while ~auto_done
      switch ncnt_auto
        case 1
          try
            %% Automatically calculate lag by interpolating depth, then time
            % STEP1 | try to interpolate depth by time (this may not work if there
            % is a time lag)
            % interpolate by time - if this is the first step, assumes no time lag
            [~,iu] = unique(ctd_date);
            lisst_interp_d = interp1(ctd_date(iu),smo_ctd(iu),lisst_time);
            % plot(gca,lisst_time,lisst_interp_d,'*')
            [~,imaxlisst] = nanmax(smo_lisst);
            [~,iminlisst] = nanmin(smo_lisst(1:imaxlisst));
            iminlisst = iminlisst+fix(imaxlisst/30);% don't choose the absolute minimum, because that could be during a soak with noisy data, so select a few meteres below
            
            dlag_depth1 = smo_lisst(iminlisst);
            dlag_depth2 = smo_lisst(imaxlisst);
            dlag1 = lisst_interp_d(iminlisst) - smo_lisst(iminlisst);
            dlag2 = lisst_interp_d(imaxlisst) - smo_lisst(imaxlisst);
            
            % Correct depth for lag assuming linear drift
            depth_bin_estimate = 1:1:fix(smo_lisst(imaxlisst));
            depth_lag_estimate = interp1([dlag_depth1 dlag_depth2],[dlag1 dlag2],depth_bin_estimate);
            depth_lag = interp1(depth_bin_estimate,depth_lag_estimate,smo_lisst);
            depth_lag = fillmissing(depth_lag,'linear','EndValues','extrap');
            depth_cor =  smo_lisst + depth_lag;
            
            %plot(gca,lisst_time,depth_cor,'*')
            
            % STEP2 | Interpolate by time by depth
            [~,imaxctd]  = max(smo_ctd);
            [~,iminctd]  = min(smo_ctd(1:imaxctd));
            
            idx_down_ctd = iminctd:imaxctd;
            [~,iuctd]    = unique(smo_ctd(idx_down_ctd));
            iuctd = idx_down_ctd(iuctd);
            
            [~,imaxlisst] = max(depth_cor);
            [~,iminlisst] = min(depth_cor(1:imaxlisst));
            
            idx_down_lisst = iminlisst:imaxlisst;
            [~,iulisst,iorig] = unique(depth_cor(idx_down_lisst));
            iulisst = idx_down_lisst(iulisst);
            
            lisst_interp_t = interp1(smo_ctd(iuctd),ctd_date(iuctd),depth_cor(iulisst));
            
            % Calculate time difference
            time_difference  = nanmedian(lisst_interp_t - lisst_time(iulisst));
            lisst_time_cor   = lisst_time + time_difference;
            time_lag_seconds = time_difference.*60.*60.*24; % convert from datenum to seconds
            
            % STEP3 | Interpolate depth by time now
            lisst_interp_d = interp1(ctd_date,smo_ctd,lisst_time_cor);
            
            [~,imaxlisst] = nanmax(lisst_interp_d);
            [~,iminlisst] = nanmin(lisst_interp_d(1:imaxlisst));
            
            dlag_depth1 = lisst_interp_d(iminlisst);
            dlag_depth2 = lisst_interp_d(imaxlisst);
            dlag1 = lisst_interp_d(iminlisst) - smo_lisst(iminlisst);
            dlag2 = lisst_interp_d(imaxlisst) - smo_lisst(imaxlisst);
            % Correct depth for lag assuming linear drift
            depth_bin_estimate = smo_lisst(iminlisst):1:smo_lisst(imaxlisst);
            depth_lag_estimate = interp1([dlag_depth1 dlag_depth2],[dlag1 dlag2],depth_bin_estimate);
            depth_lag = interp1(depth_bin_estimate,depth_lag_estimate,smo_lisst);
            depth_lag = fillmissing(depth_lag,'linear','EndValues','extrap');
            depth_cor =  smo_lisst + depth_lag;
            % STEP4 | Interpolate time by depth now
            [~,imaxlisst] = max(depth_cor);
            [~,iminlisst] = min(depth_cor(1:imaxlisst));
            
            idx_down_lisst = iminlisst:imaxlisst;
            [~,iulisst,iorig] = unique(depth_cor(idx_down_lisst));
            iulisst = idx_down_lisst(iulisst);
            
            lisst_interp_t = interp1(smo_ctd(iuctd),ctd_date(iuctd),depth_cor(iulisst));
            
            % Calculate time difference
            time_difference  = nanmedian(lisst_interp_t - lisst_time(iulisst));
            lisst_time_cor   = lisst_time + time_difference;
            time_lag_seconds = time_difference.*60.*60.*24; % convert from datenum to seconds
          catch
            
          end
          %plot(gca,lisst_time,depth_cor,'y.')
          %plot(gca,lisst_time_cor,depth_cor,'y.')
        case 2
          %% Automatically detect turning points
          try
            for profiler = 1:2
              if profiler == 1
                time = lisst_time;
                depth = smo_lisst;
              elseif profiler == 2
                time = ctd_date;
                depth = smo_ctd;
              end
              try
                TP = find_turning_points(time,depth);
              catch
                keyboard
              end
              if profiler == 1 %% LISST
                tpLISST = TP;
                [~,imaxLISST] = max(smo_lisst(tpLISST));
                tpLISST = [tpLISST(imaxLISST) find(isfinite(smo_lisst),1,'last')];
              elseif profiler == 2 %% CTD
                tpCTD   = TP;
                [~,imaxCTD] = max(smo_ctd(tpCTD));
                tpCTD = [tpCTD(imaxCTD)  find(isfinite(smo_ctd),1,'last')];
              end
              if isempty(TP)
                skip_auto = 1;
                select_manually = 1;
              end
            end
            
            
            % Calculate lag using automatically selected points
            % time lag
            time_difference  = ctd_date(tpCTD(1)) - lisst_time(tpLISST(1));
            lisst_time_cor   = lisst_time + time_difference;
            time_lag_seconds = time_difference.*60.*60.*24; % convert from datenum to seconds
            % depth lag
            dlag1       = smo_ctd(tpCTD(1)) - smo_lisst(tpLISST(1));
            dlag2       = smo_ctd(tpCTD(2)) - smo_lisst(tpLISST(2));
            dlag_depth1 = smo_lisst(tpLISST(1));
            dlag_depth2 = smo_lisst(tpLISST(2));
            
            % Correct depth for lag assuming linear drift
            depth_bin_estimate = 0:1:fix(max(smo_lisst));
            depth_lag_estimate = interp1([dlag_depth1 dlag_depth2],[dlag1 dlag2],depth_bin_estimate);
            depth_lag = interp1(depth_bin_estimate,depth_lag_estimate,smo_lisst);
            depth_lag = fillmissing(depth_lag,'linear','EndValues','extrap');
            depth_cor =  smo_lisst + depth_lag;
            
            h1 = plot(ax1,lisst_time(tpLISST), smo_lisst(tpLISST), 'kp', 'markerfacecolor', [0.4 0.4 0.4],'markersize', 10,'DisplayName','matching LISST');
            hc = plot(ax1,ctd_date(tpCTD), smo_ctd(tpCTD), 'ks', 'markerfacecolor', 'b','markersize', 10,'DisplayName','matching CTD');
          catch
            select_manually = 1;
            break
          end
        otherwise
          fprintf('NOT SET UP YET\n')
          keyboard
      end
      %% plot automatically detected points for matching ctd&lisst profiles
      h1_cor = plot(ax1,lisst_time_cor,depth_cor,'g.-','LineWidth',1,'DisplayName','LISST corrected');
      hlagt = text(ax1,0.02,0.15,['time offset  = ' num2str(time_lag_seconds,'%.2f') ' seconds'],'units','normalized','FontWeight','bold');
      hlag1 = text(ax1,0.02,0.10, ['depth offset = ' num2str(dlag1,'%.2f')  'm @ ' num2str(dlag_depth1,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
      hlag2 = text(ax1,0.02,0.05, ['depth offset = ' num2str(dlag2,'%.2f') 'm @ ' num2str(dlag_depth2,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
      
      if all(isnan(depth_cor)) == 1 
        if ncnt_auto == 1
          auto_choice_done = 1;
          ncnt_auto = ncnt_auto + 1;
       
        end
      else
        auto_choice_done = 0;
      end
      while ~auto_choice_done
        fprintf('\nCast %d - Does the lag corrected profile look good?\n',ncast)
        if ncnt_auto > 1
          fprintf('   <0>  No - select points manually\n')
        else
          fprintf('   <0>  No - try again\n')
        end
        fprintf('   <1>  Yes\n')
        fprintf('   <99> Stop\n')
        chc_auto = input('   Enter choice: ');
        if isempty(chc_auto) % Default to YES
          chc_auto = 1;
        end
        if chc_auto == 0
          if ncnt_auto == auto_methods
            select_manually  = 1;
            auto_choice_done = 1;
            auto_done = 1;
          elseif ncnt_auto == 1
            auto_choice_done = 1;
            ncnt_auto = ncnt_auto + 1;
          end
          try delete(h1); delete(hc); end
          delete(h1_cor);
          delete(hlagt);
          delete(hlag1);
          delete(hlag2);
        elseif chc_auto == 99
          keyboard
          auto_choice_done = 0;
        elseif chc_auto == 1 
          auto_done = 1;
          auto_choice_done = 1;
          select_manually = 0;
          % save figure
          if cfg.savefig
            figname = fullfile(cfg.path.dir_figs,[cfg.project '_ctd' scast '_sensor_lag']);
            standard_printfig_lowrespng(figname)
          end
          % Clear figure
          try delete(h1); delete(hc);end
          delete(h1_cor);
          delete(hlagt);
          delete(hlag1);
          delete(hlag2);
        end
      end % WHILE ~auto_choice_done
    end % WHILE ~auto_done
    
    %% Select manually - first try moving the LISST profile
    if select_manually
      try
        hlisst_cor = copyobj(hlisst,ax1);
        hlisst_cor.Color = 'g';
        hlisst_cor.LineWidth = 5;
        % move the green line a bit so it is obvious which one to select
        yshift = fix(ax1.YLim(2)/10);
        hlisst_cor.YData = hlisst_cor.YData + yshift;
        ht = text(ax1,0.02,0.95,{'Drag LISST profile (green) to corrected position'; 'Click outside axe box when finished'},...
          'Color','g','units','normalized','FontWeight','bold','FontSize',26);
        % Jean-Pascal Rueff (2021). Click-n'-drag plot (https://www.mathworks.com/matlabcentral/fileexchange/5751-click-n-drag-plot), MATLAB Central File Exchange. Retrieved June 14, 2021.
        done_moving = 0;
        while ~done_moving
          handles = interactive_move;
          done_moving = input('  Finished moving to correct position? <1/0> ');
        end
        % make sure interactive mode is stopped
        try uirestore(handles.init_state);catch; end

        %% calculate depth lag
        [~,lisst_point_deep]    = nanmax(hlisst.YData);
        [~,lisst_point_shallow] = nanmin(hlisst.YData);
        
        dlag_depth1 = hlisst.YData(lisst_point_deep);
        dlag1 = hlisst_cor.YData(lisst_point_deep) - hlisst.YData(lisst_point_deep);
        dlag_depth2 = hlisst.YData(lisst_point_shallow);
        dlag2 = hlisst_cor.YData(lisst_point_shallow) - hlisst.YData(lisst_point_shallow);
        if isequal(dlag1, dlag2)
          dlag1 = dlag1+0.10;
        end
        depth_cor = hlisst_cor.YData;
        depth_lag =  smo_lisst - depth_cor';
        %% calculate time lag
        time_difference = hlisst_cor.XData(1) - hlisst.XData(1);
        lisst_time_cor  = lisst_time + time_difference;
        time_lag_seconds = time_difference.*60.*60.*24; % convert from datenum to seconds
        
        hlagt = text(ax1,0.02,0.15,['time offset  = ' num2str(time_lag_seconds,'%.2f') ' seconds'],'units','normalized','FontWeight','bold');
        hlag1 = text(ax1,0.02,0.10, ['depth offset = ' num2str(dlag1,'%.2f')  'm @ ' num2str(dlag_depth1,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
        hlag2 = text(ax1,0.02,0.05, ['depth offset = ' num2str(dlag2,'%.2f') 'm @ ' num2str(dlag_depth2,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
        
        %% Choose whether corrected profile is okay, or if not, try again
        chc = input('  Does the corrected profile look good? <1/0> ');
        if chc == 0
          done_moving = 0;
          select_manually = 1;
        else
          done_moving = 1;
          select_manually = 0;
        end
        try
          uirestore(handles.init_state);
        catch
        end
        if done_moving
          % save figure
          if cfg.savefig
            figname = fullfile(cfg.path.dir_figs,[cfg.project '_cast' scast '_sensor_lag']);
            standard_printfig_lowrespng(figname)
          end
        end
        
        % clear figure
        delete(hlisst_cor); delete(ht);
        delete(hlag1); delete(hlag2); delete(hlagt);
      catch
        % Select manually by clicking on points
        keyboard
        select_manually = 1;
      end
    end
  
    
    %% Select corresponding deepest and shallowest points in CTD&LISST profiles to calculate the time and depth lags
    % depth lag is not constant but seems to be linear with depth so
    % calculated as such... time lag is constant
    if select_manually
      done_selecting_points = 0; % select points on graph
      while ~done_selecting_points
        %% Find deepest depth offset and time offset
        fprintf('-------------------------------------\n')
        fprintf('-------------------------------------\n')
        fprintf('  Select DEEPEST CTD profile point\n')
        ht = text(0.2,0.95,'Select DEEPEST CTD profile point','FontWeight','bold','FontSize',30,'Units','normalized','Color','r');
        ctd_pt = ginput(1);
        ctd_point_deep = find(ctd_date >= ctd_pt(1),1);
        hp1 = plot( ctd_date(ctd_point_deep), smo_ctd(ctd_point_deep), 'kp', 'markerfacecolor', 'g','markersize',15);
        delete(ht);
        
        fprintf('  Select corresponding DEEPEST LISST profile point\n')
        ht = text(0.2,0.95,'Select corresponding DEEPEST LISST profile point','FontWeight','bold','FontSize',30,'Units','normalized','Color','r');
        lisst_pt = ginput(1);
        lisst_point_deep = find(lisst_time >= lisst_pt(1),1);
        hp2 = plot( lisst_time(lisst_point_deep), smo_lisst(lisst_point_deep), 'kp', 'markerfacecolor', 'r','markersize',15);
        delete(ht);
        
        %% Find shallowest depth offset
        fprintf('  Select SHALLOW CTD profile point\n')
        ht = text(0.2,0.95,'Select SHALLOW CTD profile point','FontWeight','bold','FontSize',30,'Units','normalized','Color','r');
        ctd_pt = ginput(1);
        ctd_point_shallow = find(ctd_date >= ctd_pt(1),1);
        if  isempty(ctd_point_shallow)
          ctd_point_shallow = find(isfinite(smo_ctd),1,'last');
        elseif isnan(smo_ctd(ctd_point_shallow))
          ctd_point_shallow = find(isfinite(smo_ctd),1,'last');
        end
        
        
        hp3 = plot( ctd_date(ctd_point_shallow), smo_ctd(ctd_point_shallow), 'kp', 'markerfacecolor', 'g','markersize',15,'DisplayName','matching CTD');
        delete(ht);
        
        fprintf('  Select corresponding SHALLOW LISST profile point\n')
        ht = text(0.2,0.95,'Select corresponding SHALLOW LISST profile point','FontWeight','bold','FontSize',30,'Units','normalized','Color','r');
        lisst_pt = ginput(1);
        lisst_point_shallow = find(lisst_time >= lisst_pt(1),1);
        if  isempty(lisst_point_shallow)
          lisst_point_shallow = find(isfinite(smo_lisst),1,'last');
        elseif isnan(smo_lisst(lisst_point_shallow))
          lisst_point_shallow = find(isfinite(smo_lisst),1,'last');
        end
        hp4 = plot( lisst_time(lisst_point_shallow), smo_lisst(lisst_point_shallow), 'kp', 'markerfacecolor', 'r','markersize',15,'DisplayName','matching LISST');
        delete(ht);
        
        %% calculate depth lag
        dlag_depth1 = smo_lisst(lisst_point_deep);
        dlag1 = smo_ctd(ctd_point_deep) - smo_lisst(lisst_point_deep);
        dlag_depth2 = smo_lisst(lisst_point_shallow);
        dlag2 = smo_ctd(ctd_point_shallow) - smo_lisst(lisst_point_shallow);
        
        %% calculate time lag
        time_difference = ctd_date(ctd_point_deep) - lisst_time(lisst_point_deep);
        lisst_time_cor = lisst_time + time_difference;
        time_lag_seconds = time_difference.*60.*60.*24; % convert from datenum to seconds
        
        % display lag estimates on the plot
        hlagt = text(ax1,0.02,0.15,['time offset  = ' num2str(time_lag_seconds,'%.2f') ' seconds'],'units','normalized','FontWeight','bold');
        hlag1 = text(ax1,0.02,0.10,['depth offset = ' num2str(dlag1,'%.2f')  'm @ ' num2str(dlag_depth1,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
        hlag2 = text(ax1,0.02,0.05,['depth offset = ' num2str(dlag2,'%.2f')  'm @ ' num2str(dlag_depth2,'%.2f') 'm'], 'units','normalized','FontWeight','bold');
        
        %% Calculate depth lag assuming linear drift & Plot correct profile to check if it looks good
        try
          depth_bin_estimate = 1:1:fix(max(smo_lisst));
          depth_lag_estimate = interp1([dlag_depth1 dlag_depth2],[dlag1 dlag2],depth_bin_estimate);
          depth_lag_estimate = fillmissing(depth_lag_estimate,'linear','EndValues','extrap');
          depth_lag = interp1(depth_bin_estimate,depth_lag_estimate,smo_lisst);
          depth_lag = fillmissing(depth_lag,'linear','EndValues','extrap');
          depth_cor =  smo_lisst + depth_lag;
        catch
          fprintf('could not calculate depth lag and corrected depth\n')
          keyboard
        end
        h1_cor = plot(ax1,lisst_time_cor,depth_cor,'g.-','LineWidth',1,'DisplayName','LISST corrected');
        %% Choose whether corrected profile is okay, or if not, try again
        chc = input('  Does the corrected profile look good? <1/0> ');
        if chc == 0
          done_selecting_points = 0;
        else
          done_selecting_points = 1;
        end
        if done_selecting_points
          % save figure
          if cfg.savefig
            figname = fullfile(cfg.path.dir_figs,[cfg.project '_cast' scast '_sensor_lag']);
            standard_printfig_lowrespng(figname)
          end
        end
        % clear figure
        delete(h1_cor)
        delete(hp1);   delete(hp2);   delete(hp3);  delete(hp4);
        delete(hlag1); delete(hlag2); delete(hlagt);
        
      end % DONE_SELECTING_POINTS
    end %% SELECT MANUALLY
    
    % Write to file: cast,datfile,time_lag,depth1,depth2,depth1_lag,depth2_lag
    fprintf(fileID,'\n%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f\n',ncast,datname,time_lag_seconds,dlag_depth1,dlag_depth2,dlag1,dlag2);
    
  end %% LOADED LAGS OR CALCULATED LAGS
  % store values to structure
  cfg.proc_options.time_lag(nf)   = time_lag_seconds;
  cfg.proc_options.depth1(nf)     = dlag_depth1;
  cfg.proc_options.depth2(nf)     = dlag_depth2;
  cfg.proc_options.depth1_lag(nf) = dlag1;
  cfg.proc_options.depth2_lag(nf) = dlag2;
  
  % store corrected time and date vectors to data_pre structure
  data_pre.depth(idx_lisst) = data_pre.depth(idx_lisst) + depth_lag;
  data_pre.date(idx_lisst)  = lisst_time_cor;
  % clear and reset axes
  cla; ax1 = gca;
end %% LOOP THROUGH CASTS
fprintf('DONE finding time and depth lags\n')

%% Diagnostic plots
% Pull out time at the beginning of each cast
times = nan(numel(cfg.proc_options.datfile),1);
for nc = 1:numel(times)

  idx_lisst = find(data_pre.cast == cfg.proc_options.cast(nc));
  times(nc) = data_pre.date(idx_lisst(1));
end
close; ax = makefig_subplots(1,2); ax1 = ax(2); % original

plot(ax1,times,cfg.proc_options.time_lag./60.0,'*-','MarkerSize',10,'LineWidth',2)
ax1.YLabel.String = 'Time lag [minutes]';
datetick(ax1,'x','keepticks','keeplimits')
grid(ax1,'on'); ax1.XTickLabel = [];
% plot depth lag at 20m and at 100m
lag_20m  = nan(numel(times),1);
lag_100m = nan(numel(times),1);
for nc = 1:numel(times)
  idx_lisst = find(data_pre.cast == cfg.proc_options.cast(nc));
  w20  = find(data_pre.depth(idx_lisst) >= 20,1);
  w100 = find(data_pre.depth(idx_lisst) >= 100,1);
  if ~isempty(w20);  lag_20m(nc)  = data_pre.depth(idx_lisst(w20))  - data_pre.depth_orig(idx_lisst(w20)); end
  if ~isempty(w100); lag_100m(nc) = data_pre.depth(idx_lisst(w100)) - data_pre.depth_orig(idx_lisst(w100)); end
end
ax2 = ax(1);
plot(ax2,times,lag_20m,'rd-','MarkerSize',6,'LineWidth',2,'DisplayName','offset @ 20m')
hold(ax2,'on'); grid(ax2,'on');
plot(ax2,times,lag_100m,'ks-','MarkerSize',6,'LineWidth',2,'DisplayName','offset @ 100m')
ax2.XLim = ax1.XLim;
datetick(ax2,'x','keepticks','keeplimits')
ax2.XLabel.String = 'date';
ax2.YLabel.String = 'Depth lag [m]';
legend(ax2,'show','location','nw')
title({strrep(cfg.project,'_','\_'); 'Time and depth offsets compared to CTD'})
% save figure
if cfg.savefig
  figname = fullfile(cfg.path.dir_figs,[cfg.project '_time_depth_offsets']);
  standard_printfig_lowrespng(figname)
end


end %% MAIN FUNCTION
