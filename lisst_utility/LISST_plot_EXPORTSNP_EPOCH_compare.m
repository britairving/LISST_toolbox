function LISST_plot_EXPORTS_EPOCH_compare(RR,SR)

% 
%%
%% Spherical average filtered seawater background scatter files 
ftypes = {'RogerRevelle' 'SallyRide'}; % file type 
ptype  = 'spherical medianMilliQ'; % processing type
if nargin < 2
  % RR = load('F:\LISST_Data\LISST_sn4041_2018_exports_np_rr1813\LISST_sn4041_2018_exports_np_rr1813_avgFiltSWzscat_spherical_qc_data_grid_ignoreflagged.mat');
  % SR = load('F:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\LISST_sn4025_2018_exports_np_sr1812_avgFSWzscat_spherical_qc_data_grid_ignoreflagged.mat');
  RR = load('/Volumes/toshiba1/LISST_Data_Structured/LISST_sn4041_2018_EXPORTS_RR201813/proc/LISST_sn4041_2018_EXPORTS_RR201813_proc_spherical_gridded.mat');
  SR = load('/Volumes/toshiba1/LISST_Data_Structured/LISST_sn4025_2018_EXPORTS_SR201812/proc/LISST_sn4025_2018_EXPORTS_SR201812_proc_spherical_gridded.mat');
end
%% 0 | Set up script options
depth_var = 'depth';% 'depth_orig'; % 'depth'
date_var  = 'date'; % 'date_orig';  % 'date'
script_opt.axcolor = [0.8 0.80 0.8];
script_opt.savefig = 1;
script_opt.savedir = '/Volumes/toshiba1/LISST_Data_Structured/LISST_sn4025_2018_EXPORTS_SR201812';%'F:\LISST_Data\';
% casts = [10, 68];
SRcasts = [20 87 130];
RRcasts = [9 52 79];

%% Loop through epochs and plot casts
for epoch = 1:3
  %% Plot total volume concentration, mean size, and surface area
  ax = makefig_subplots(3,1); set(gcf,'Name','Exports Epoch Comparison Cast');
  
  for ndata = 1:2
    if ndata == 1
      opt    = RR.cfg;
      data_grid = RR.data_grid;
      wcast  = find(data_grid.cast == RRcasts(epoch));
      RRwcast = wcast;
      clr = 'k';
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
    elseif ndata == 2
      opt    = SR.cfg;
      data_grid = SR.data_grid;
      wcast  = find(data_grid.cast == SRcasts(epoch));
      SRwcast = wcast;
      clr = 'b';
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
    end
    CastID = ['Cast' num2str(data_grid.cast(wcast))];
    lstring = ftypes{ndata};
    % Load BRU & DAT files
    fprintf(' Plotting %s LISST %s data\n',lstring, CastID)
    dias    = opt.proc_options.dias;     % opt.meta.dias = 'midpoint of size bins'
    dep     = data_grid.(depth_var);
    if numel(zscfile) > 20
      zscfile = erase(zscfile,[opt.project '_']);
    end
    lstring = [strrep(lstring,'_','\_') ' ' strrep(zscfile,'_','\_')];
    %% Plot profile statistics
    %max_depth = max(data_grid.(depth_var));
    ylims = [0 200]; %max_depth];

    x = data_grid.tot_vol_concentration(wcast,:);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(1),x,y,ylims,clr,lstring);
    ax(1).XLabel.String = 'Total Volume Concentration [µl/l]';
    hold(ax(1),'on');
    
    x = data_grid.mean_size(wcast,:);
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(2),x,y,ylims,clr,lstring);
    ax(2).XLabel.String = 'Mean Size [µm]';
    ax(2).YLabel = []; ax(2).YTickLabel = [];
    hold(ax(2),'on');
    
    x = data_grid.surf_area(wcast,:); % squeeze(data_grid.cscat(wcast,:,:));
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(3),x,y,ylims,clr,lstring);
    ax(3).XLabel.String = 'Surface Area [cm^2/l]';
    ax(3).YLabel = []; ax(3).YTickLabel = [];  ax(3).Title.String = [];
    hold(ax(3),'on');
  end
  ax(2).Title.String = {['Sally Ride Cast ' num2str(SRcasts(epoch)) ' [' num2str(SR.data_grid.lat(SRwcast),'%.2f') '\circN, ' num2str(SR.data_grid.lon(SRwcast),'%.2f') '\circW]'];...
                       ['Roger Revelle Cast ' num2str(RRcasts(epoch)) ' [' num2str(RR.data_grid.lat(RRwcast),'%.2f') '\circN, ' num2str(RR.data_grid.lon(RRwcast),'%.2f') '\circW]']};
  legend(ax(1),'show','location','se')
  legend(ax(2),'show','location','sw')
  legend(ax(3),'show','location','se')
  
  if script_opt.savefig
    figname = fullfile(script_opt.savedir,['LISST_Exports_compare_' ptype '_Epoch' num2str(epoch)]);
    standard_printfig_highrespng(figname)
  else
    keyboard
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
end