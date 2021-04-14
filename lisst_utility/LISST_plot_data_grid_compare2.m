% dat1 = load('X:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\march10\LISST_sn4025_2018_exports_np_sr1812_binned_L2270456_average300_to_1000m.mat');
% dat2 = load('X:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\march10\LISST_sn4025_2018_exports_np_sr1812_binned_zscat_MilliQ_avg.mat');
% dat3 = load('X:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\march10\LISST_sn4025_2018_exports_np_sr1812_binned_zscat_FSW_avg.mat');

%% RR1813 Revelle
% ftypes = {'nonspherical' 'spherical'};
% dat1 = load('F:\LISST_Data\LISST_sn4041_2018_exports_np_rr1813\LISST_sn4041_2018_exports_np_rr1813_nearestFiltSWzscat_nonspherical_qc_binned_ignoreflagged.mat');
% dat2 = load('F:\LISST_Data\LISST_sn4041_2018_exports_np_rr1813\LISST_sn4041_2018_exports_np_rr1813_nearestFiltSWzscat_spherical_qc_binned_ignoreflagged.mat');
% SR1812 Sally Ride
% ftypes = {'nonspherical' 'spherical'};
% dat1 = load('F:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\LISST_sn4025_2018_exports_np_sr1812_nearestFSWzscat_nonspherical_qc_binned_ignoreflagged.mat');
% dat2 = load('F:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\LISST_sn4025_2018_exports_np_sr1812_nearestFSWzscat_spherical_qc_binned_ignoreflagged.mat');

%% SR1812 backgrounds
% ftypes = {'factory_zscat' 'insitu_zscat' 'median_FSW_zscat' 'median_MilliQ_zscat'};
% dat1 = load('/Volumes/toshiba1/LISST_Data/LISST_sn4025_2018_exports_np_sr1812/LISST_sn4025_2018_exports_np_sr1812_spherical_factory_zscat_binned.mat');
% dat2 = load('/Volumes/toshiba1/LISST_Data/LISST_sn4025_2018_exports_np_sr1812/LISST_sn4025_2018_exports_np_sr1812_spherical_insitu_zscat_binned.mat');
% dat3 = load('/Volumes/toshiba1/LISST_Data/LISST_sn4025_2018_exports_np_sr1812/LISST_sn4025_2018_exports_np_sr1812_spherical_median_FSW_zscat_binned.mat');
% dat4 = load('/Volumes/toshiba1/LISST_Data/LISST_sn4025_2018_exports_np_sr1812/LISST_sn4025_2018_exports_np_sr1812_spherical_median_MilliQ_zscat_binned.mat');



%% 0 | Set up script options
cfg = dat1.cfg; % all the same - just use first
depth_var = 'depth';% 'depth_orig'; % 'depth'
date_var  = 'date'; % 'date_orig';  % 'date'
script_opt.axcolor = [0.8 0.80 0.8];
script_opt.savefig = 1;
casts = dat1.data_grid.cast;
% casts = [10, 68];

num_casts = numel(casts);
num_types = numel(ftypes);
clrs = jet(num_types);

plot_tot_PSD  = 1;
plot_avg_cast = 1;
plot_all_cast = 0;
%% scattering 
% scat(i,:)=data(i,1:32)/tau(i)-(zsc(1:32)*data(i,36)/zsc(36));

if plot_avg_cast
  fprintf(' Plotting LISST %s data\n',cfg.project)
  
  %% Plot total volume concentration, mean size, and surface area, and total particle concetration
  if plot_tot_PSD
    ax = makefig_subplots(4,1); set(gcf,'Name',dat1.cfg.project);
  else
    ax = makefig_subplots(3,1); set(gcf,'Name',dat1.cfg.project);
  end
  
  for ndata = 1:num_types
    if ndata == 1
      cfg     = dat1.cfg;
      dat_bin = dat1.data_grid;
      [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
      [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
    elseif ndata == 2
      cfg     = dat2.cfg;
      dat_bin = dat2.data_grid;
      [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
      [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
    elseif ndata == 3
      cfg     = dat3.cfg;
      dat_bin = dat3.data_grid;
      [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
      [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
    elseif ndata == 4
      cfg     = dat4.cfg;
      dat_bin = dat4.data_grid;
      [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
      [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
    end
    
    lstring = strrep(ftypes{ndata},'_','\_');
    dias    = cfg.dias;     % cfg.meta.dias = 'midpoint of size bins'
    dep     = dat_bin.(depth_var);
    
    %% Plot profile statistics
    %max_depth = max(dat_bin.(depth_var));
    ylims = [0 200]; %max_depth];
    
    % total volume concentration
    
    x = nanmean(dat_bin.tot_vol_concentration,1);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(1),x,y,ylims,clrs(ndata,:),lstring);
    ax(1).XLabel.String = 'Total Volume Concentration [µl/l]';
    hold(ax(1),'on');
    
    % Mean Size
    x = nanmean(dat_bin.mean_size,1);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(2),x,y,ylims,clrs(ndata,:),lstring);
    ax(2).XLabel.String = 'Mean Size [µm]';
    ax(2).YLabel = []; ax(2).YTickLabel = [];
    hold(ax(2),'on');
    
    % Surface Area
    x = nanmean(dat_bin.surf_area,1);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(3),x,y,ylims,clrs(ndata,:),lstring);
    ax(3).XLabel.String = 'Surface Area [cm^2/l]';
    ax(3).YLabel = []; ax(3).YTickLabel = [];  ax(3).Title.String = [];
    hold(ax(3),'on');
    
    % total particle concentration
    if plot_tot_PSD
      % calculate total particle size distribution
      x = nansum(dat_bin.PSD,3);
      x = nanmean(x,1);
      igood = isfinite(x); x = x(igood); y = dep(igood);
      plot_against_pressure(ax(4),x,y,ylims,clrs(ndata,:),lstring);
      ax(4).XLabel.String = 'Total Particle Concentration [#/µm^3/m]';
      ax(4).YLabel = []; ax(4).YTickLabel = [];  ax(4).Title.String = [];
      hold(ax(4),'on');
    end
    
    % Add title
    ax(1).Title.String = {[strrep(cfg.project,'_','\_') ' | Averaged casts ' num2str(casts(1)) '-' num2str(casts(end))]};
    ax(1).Title.FontSize = 20;
    ax(1).Title.HorizontalAlignment = 'left';
  end
  legend(ax(1),'show','location','se')
  legend(ax(2),'show','location','se')
  legend(ax(3),'show','location','se')
  if plot_tot_PSD
    legend(ax(4),'show','location','se')
  end
  if script_opt.savefig
    figname = fullfile(cfg.path.figs,[cfg.project '_averagedcasts' num2str(casts(1)) 'to' num2str(casts(end)) '_compare']);
    standard_printfig_highrespng(figname)
  else
    keyboard
  end
end % Plot average cast

if plot_all_cast
  for nc = 1:num_casts
    wcast = dat1.data_grid.cast == casts(nc);
    CastID = ['Cast' num2str(dat1.data_grid.cast(wcast))];
    % Load BRU & DAT files
    fprintf(' Plotting LISST %s data\n',CastID)
    
    %% Plot total volume concentration, mean size, and surface area
    if plot_tot_PSD
      ax = makefig_subplots(4,1); set(gcf,'Name',dat1.cfg.project);
    else
      ax = makefig_subplots(3,1); set(gcf,'Name',dat1.cfg.project);
    end
    
    for ndata = 1:num_types
      if ndata == 1
        cfg     = dat1.cfg;
        dat_bin = dat1.data_grid;
        [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
        [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
      elseif ndata == 2
        cfg     = dat2.cfg;
        dat_bin = dat2.data_grid;
        [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
        [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
      elseif ndata == 3
        cfg     = dat3.cfg;
        dat_bin = dat3.data_grid;
        [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
        [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
      elseif ndata == 4
        cfg     = dat4.cfg;
        dat_bin = dat4.data_grid;
        [~,datfile,~] = fileparts(dat_bin.datfile{wcast});
        [~,zscfile,~] = fileparts(dat_bin.zscfile{wcast});
      end
      
      lstring = strrep(ftypes{ndata},'_','\_');
      dias    = cfg.dias;     % cfg.meta.dias = 'midpoint of size bins'
      dep     = dat_bin.(depth_var);
      
      %% Plot profile statistics
      %max_depth = max(dat_bin.(depth_var));
      ylims = [0 200]; %max_depth];
      
      % total volume concentration
      x = dat_bin.tot_vol_concentration(wcast,:);
      igood = isfinite(x); x = x(igood); y = dep(igood);
      plot_against_pressure(ax(1),x,y,ylims,clrs(ndata,:),lstring);
      ax(1).XLabel.String = 'Total Volume Concentration [µl/l]';
      hold(ax(1),'on');
      
      % Mean Size
      x = dat_bin.mean_size(wcast,:);
      igood = isfinite(x); x = x(igood); y = dep(igood);
      plot_against_pressure(ax(2),x,y,ylims,clrs(ndata,:),lstring);
      ax(2).XLabel.String = 'Mean Size [µm]';
      ax(2).YLabel = []; ax(2).YTickLabel = [];
      hold(ax(2),'on');
      
      % Surface Area
      x = dat_bin.surf_area(wcast,:); % squeeze(dat_bin.cscat(wcast,:,:));
      igood = isfinite(x); x = x(igood); y = dep(igood);
      plot_against_pressure(ax(3),x,y,ylims,clrs(ndata,:),lstring);
      ax(3).XLabel.String = 'Surface Area [cm^2/l]';
      ax(3).YLabel = []; ax(3).YTickLabel = [];  ax(3).Title.String = [];
      hold(ax(3),'on');
      
      % total particle concentration
      if plot_tot_PSD
        % calculate total particle size distribution
        x = nansum(dat_bin.PSD(wcast,:,:),3);
        igood = isfinite(x); x = x(igood); y = dep(igood);
        plot_against_pressure(ax(4),x,y,ylims,clrs(ndata,:),lstring);
        ax(4).XLabel.String = 'Total Particle Concentration [#/µm^3/m]';
        ax(4).YLabel = []; ax(4).YTickLabel = [];  ax(4).Title.String = [];
        hold(ax(4),'on');
      end
      
      % Add title
      ax(1).Title.String = {[strrep(cfg.project,'_','\_') ' ' strrep(CastID,'_', ' ') ': ' datfile]};
      ax(1).Title.FontSize = 20;
      ax(1).Title.HorizontalAlignment = 'left';
      
    end
    legend(ax(1),'show','location','se')
    legend(ax(2),'show','location','se')
    legend(ax(3),'show','location','se')
    if plot_tot_PSD
      legend(ax(4),'show','location','se')
    end
    if script_opt.savefig
      figname = fullfile(cfg.path.figs,[cfg.project '_' CastID '_compare']);
      standard_printfig_highrespng(figname)
    else
      keyboard
    end
  end
  
end


%   semilogy(x,y,clr)
%% FUNCTION plot_against_pressure
function plot_against_pressure(ax,x,y,ylim,clr,legend_string)
if nargin < 5; clr = 'k'; end % default to black
axes(ax); hold(ax,'on'); grid(ax,'on')% make sure proper axes is current
if nargin == 6
  plot(ax,x,y,'-','color',clr,'LineWidth',2,'DisplayName',legend_string);
else
   plot(ax,x,y,'-','color',clr,'LineWidth',2)
end
set(ax,'XAxisLocation','bottom','YDir','reverse');
ylabel(ax,'Depth [m]')
set(ax,'YLim',ylim)
set(ax,'Color',[0.8 0.8 0.8])
end %% FUNCTION plot_against_pressure
