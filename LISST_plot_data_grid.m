function LISST_plot_data_grid(cfg,data_grid)
%FUNCTION LISST_plot_data_grid
%
%  Syntax:
%    LISST_plot_data_grid(cfg,data_grid,meta_proc)
%  
%  Description:
%    
%
%  Refereces:
%    % http://www.sequoiasci.com/article/interpreting-particle-data/
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
fprintf('This could use some attention to make it more useful...\n')

%% 0 | Set up script options
zname = unique(data_grid.zscfile);
if numel(zname) == 1
  zname = [char(erase(zname,'.asc')) '_'];
else
  zname = ''; % multiple background files used
end
script_opt.axcolor = [0.8 0.80 0.8];
script_opt.plot_PSD  = 0;
script_opt.plot_scat = 0;
script_opt.plot_scat2 = 0;
script_opt.plot_3vars = 1; 
ucasts = unique(data_grid.cast);
% ucasts = [10, 68];

num_casts = numel(ucasts);
dep       = data_grid.depth;
if script_opt.plot_PSD
  for nc = 1:num_casts
    wcast = find(data_grid.cast == ucasts(nc));
    CastID = ['Cast' num2str(data_grid.cast(wcast(1)))];
    % Load BRU & DAT files
    fprintf(' Plotting LISST %s data\n',CastID)
    [~,datfile,~] = fileparts(data_grid.datfile{wcast(1)});
    [~,zscfile,~] = fileparts(data_grid.zscfile{wcast(1)});
    %% Plot Concentration Size Distribution (CSDn) at all depths within one station
    makefig;ax = gca; hold(ax,'on');
    i=0;
    dep  = data_grid.depth;
    depths_to_show = [1:2:50 50:10:200];
    cmap = lansey(numel(depths_to_show));
    for k = depths_to_show
      %     idepth = data_grid.depth == depths_to_show(k);
      i=i+1;
      PSD = squeeze(data_grid.PSD(wcast,k,:));
      dias_mid = cfg.proc_options.dias;
      iz    = find(PSD == 0);
      PSD(iz) = [];
      dias_mid(iz)  = [];
      fprintf('%s\n',[num2str(dep(k)) 'm'])
      % plot(szav,CSDnk,'Color',cmap_profbysize(i,:),'Marker','o','LineStyle','-','DisplayName',[num2str(round(dep(k))) 'm'])
      plot(dias_mid,PSD,'Color',cmap(i,:),'Marker','o','LineStyle','-','LineWidth',2,'DisplayName',[num2str(dep(k)) 'm'])
    end
    try
      axes(ax);
      hl = legend(gca,'show'); hl.FontSize = 10;
    catch
    end
    ax.YDir = 'normal';
    set(gca,'Color',script_opt.axcolor)
    set(gca, 'YScale', 'log', 'XScale','log');%,'Xlim',[0.05,25],'YLim',[10^-3,10^4]);
    xlabel('Diameter mid point [µm]');
    ylabel('Particle Size Distribution [#/µm^3/m]');
    hold off
    if isfield(data_grid,'station')
      title({['Station: ' data_grid.station{wcast} ' | ' strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]});
    else
      title({[strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]});
    end
    
    text(gca,0.02,0.97,'PSDs by depth','units','normalized','FontSize',24)
    
    set(gcf,'Name',[strrep(cfg.project,'_',' ') ' PSDs by depth'])
    if cfg.savefig
      figname = fullfile(cfg.path.dir_figs,[cfg.project '_' zname CastID 'PSD']);
      standard_printfig_lowrespng(figname)
    else
      keyboard
    end
  end
end


if script_opt.plot_scat
  for nc = 1:num_casts
    wcast = data_grid.cast == ucasts(nc);
    CastID = ['Cast' num2str(data_grid.cast(wcast))];
    fprintf(' Plotting LISST %s data\n',CastID)
    [~,datfile,~] = fileparts(cfg.proc_options.datfile{wcast});
    [~,zscfile,~] = fileparts(cfg.proc_options.zscfile{wcast});
    %% scattering  at all depths within one station
    makefig;ax = gca; hold(ax,'on');
    dep  = data_grid.depth;
    cmap = lansey(32);
    for n = 1:32
      plot(ax, squeeze(data_grid.scat(wcast,:,n)),dep,'Color',cmap(n,:),'Marker','s','LineStyle','-','LineWidth',2,'DisplayName',[num2str(cfg.proc_options.dias(n),'%.1f') 'µm'])
      %plot(ax, squeeze(data_grid.cscat(wcast,n)),dep,'Color',cmap(n,:),'Marker','s','LineStyle','-','LineWidth',2,'DisplayName',[num2str(cfg.proc_options.dias(n),'%.1f') 'µm'])
    end
    try
      hl = legend(ax,'show'); hl.FontSize = 10;
    catch
    end
    set(gca,'Color',script_opt.axcolor)
    ax.XScale = 'linear';
    
    ax.YDir   = 'rev';
    ylabel('Depth [m]');
    xlabel('Scattering [digital counts]');
    hold off
    if isfield(data_grid,'station')
      title({['Station: ' data_grid.station{wcast} ' | ' strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]});
    else
      title({[strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]});
    end
    
    %text(gca,0.02,0.97,'scattering by depth','units','normalized','FontSize',24)

    set(gcf,'Name',[strrep(cfg.project,'_',' ') ' scat by depth'])

    if cfg.savefig
      figname = fullfile(cfg.path.dir_figs,[cfg.project '_' zname CastID 'scat']);
      standard_printfig_lowrespng(figname)
    else
      keyboard
    end
  end
end

if script_opt.plot_scat2
  for nc = 1:num_casts
    wcast = find(data_grid.cast == ucasts(nc));
    CastID = ['Cast' num2str(data_grid.cast(wcast(1)))];
    fprintf(' Plotting LISST %s data\n',CastID)
    [~,datfile,~] = fileparts(data_grid.datfile{wcast(1)});
    [~,zscfile,~] = fileparts(data_grid.zscfile{wcast(1)});
    %% scattering  at all depths within one station
    makefig;ax = gca; hold(ax,'on');
    dep  = data_grid.depth;
    cmap = lansey(32);
    plot(ax,1:32, squeeze(data_grid.scat(wcast,:,:)),'LineStyle','-','LineWidth',2)
      %plot(ax, squeeze(data_grid.cscat(wcast,n)),dep,'Color',cmap(n,:),'Marker','s','LineStyle','-','LineWidth',2,'DisplayName',[num2str(cfg.proc_options.dias(n),'%.1f') 'µm'])
    set(gca,'Color',script_opt.axcolor)
    ax.XScale = 'linear';
    
    ax.YDir   = 'rev';
    xlabel('Ring Detector no.');
    ylabel('Scattering [digital counts]');
    hold off
    
    if isfield(data_grid,'station')
      title({['Station: ' data_grid.station{wcast} ' | ' strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]});
    else
      title({[strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]});
    end
    %text(gca,0.02,0.97,'scattering by depth','units','normalized','FontSize',24)
    
    set(gcf,'Name',[strrep(cfg.project,'_',' ') ' scat by depth'])
    if cfg.savefig
      figname = fullfile(cfg.path.dir_figs,[cfg.project '_' zname CastID 'scat']);
      standard_printfig_lowrespng(figname)
    else
      keyboard
    end
  end
end


if script_opt.plot_3vars
  for nc = 1:num_casts
    wcast = find(data_grid.cast == ucasts(nc));
    CastID = ['Cast' num2str(data_grid.cast(wcast))];
    % Load BRU & DAT files
    fprintf(' Plotting LISST %s data\n',CastID)
    [~,datfile,~] = fileparts(data_grid.datfile{wcast(1)});
    [~,zscfile,~] = fileparts(data_grid.zscfile{wcast(1)});
    %% Plot profile statistics
    ylims = [0 200]; %max_depth];
    
    %% Plot total volume concentration, mean size, and surface area
    ax = makefig_subplots(3,1); set(gcf,'Name',cfg.project);
    
    x = data_grid.tot_vol_concentration(wcast,:);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(1),x,y,ylims,'k');
    ax(1).XLabel.String = 'Total Volume Concentration [µl/l]';
    
    x = data_grid.mean_size(wcast,:);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(2),x,y,ylims,'k');
    ax(2).XLabel.String = 'Mean Size [µm]';
    ax(2).YLabel = []; ax(2).YTickLabel = [];
    
    if isfield(data_grid,'station')
      ax(2).Title.String = {['Station: ' data_grid.station{wcast} ' | ' strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]};
    else
      ax(2).Title.String = {[strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]};
    end
    
    
    
    x = data_grid.surf_area(wcast,:);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(3),x,y,ylims,'k');
    ax(3).XLabel.String = 'Surface Area [cm^2/l]';
    ax(3).YLabel = []; ax(3).YTickLabel = [];  ax(3).Title.String = [];
    
    if cfg.savefig
      figname = fullfile(cfg.path.figs,[cfg.project '_' CastID '_vsdepth']);
      standard_printfig_lowrespng(figname)
    else
      keyboard
    end
  end
end
      
%% FUNCTION plot_against_pressure
  function plot_against_pressure(ax,x,y,ylim,clr)
    if nargin < 5; clr = 'k'; end % default to black
    axes(ax); hold(ax,'on'); grid(ax,'on')% make sure proper axes is current
    plot(ax,x,y,'-','color',clr,'LineWidth',2);
    set(ax,'XAxisLocation','bottom','YDir','reverse');
    ylabel(ax,'Pressure [dbar]')
    if nargin > 3
      set(ax,'YLim',ylim)
    else
      set(ax,'YLim',[0 script_opt.maxdepth])
    end
    set(ax,'Color',script_opt.axcolor)
  end %% FUNCTION plot_against_pressure


end %% MAIN FUNCTION 