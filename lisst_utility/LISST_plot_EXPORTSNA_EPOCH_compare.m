function LISST_plot_EXPORTSNA_EPOCH_compare(JC,DY)

% 
%%
%% Spherical average filtered seawater background scatter files 
ftypes = {'JamesCook' 'Discovery'}; % file type 
ptype  = 'spherical medianMQWvsFactoryZSCAT'; % processing type
if nargin < 2
  DY = load('F:\LISST_Data_Structured\LISST_sn4025_2021_EXPORTS_DY131\proc\LISST_sn4025_2021_EXPORTS_DY131_proc_spherical_gridded.mat');
  JC = load('F:\LISST_Data_Structured\LISST_sn4041_2021_EXPORTS_JC214\proc\LISST_sn4041_2021_EXPORTS_JC214_proc_spherical_gridded.mat');
  %   DY = load('C:\Users\Brita Irving\Downloads\drive-download-20220818T013420Z-001\LISST_sn4025_2021_EXPORTS_DY131_proc_spherical_gridded.mat');
  %   JC = load('C:\Users\Brita Irving\Downloads\drive-download-20220818T013420Z-001\LISST_sn4041_2021_EXPORTS_JC214_proc_spherical_gridded.mat');
end
%% 0 | Set up script options
depth_var = 'depth';% 'depth_orig'; % 'depth'
date_var  = 'date'; % 'date_orig';  % 'date'
script_opt.axcolor = [0.8 0.80 0.8];
script_opt.savefig = 1;
script_opt.savedir = 'F:\LISST_Data_Structured\';%'F:\LISST_Data\';

DYcasts = [19 31 82 94]; % Disovery found on "Discovery Matchups.xlsx"
JCcasts = [12 23 53 57]; % JamesCook found on "Cook Matchups.xlsx"
%% Loop through epochs and plot casts
for epoch = 1:numel(DYcasts)
  %% Plot total volume concentration, mean size, and surface area
  ax = makefig_subplots(3,1); set(gcf,'Name','Exports North Atlantic Epoch Comparison Cast');
  
  for ndata = 1:2
    if ndata == 1
      opt    = JC.cfg;
      data_grid = JC.data_grid;
      wcast  = find(data_grid.cast == JCcasts(epoch));
      JCwcast = wcast;
      clr = 'k';
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
    elseif ndata == 2
      opt    = DY.cfg;
      data_grid = DY.data_grid;
      wcast  = find(data_grid.cast == DYcasts(epoch));
      DYwcast = wcast;
      clr = 'b';
      try
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
      catch
        keyboard
      end
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
  ax(2).Title.String = {['Discovery Cast ' num2str(DYcasts(epoch)) ' [' num2str(DY.data_grid.lat(DYwcast),'%.2f') '\circN, ' num2str(DY.data_grid.lon(DYwcast),'%.2f') '\circW]'];...
                       ['James Cook Cast ' num2str(JCcasts(epoch)) ' [' num2str(JC.data_grid.lat(JCwcast),'%.2f') '\circN, ' num2str(JC.data_grid.lon(JCwcast),'%.2f') '\circW]']};
  hl1 = legend(ax(1),'show','location','se');
  hl2 = legend(ax(2),'show','location','sw');
  hl3 = legend(ax(3),'show','location','se');
  hl1.FontSize = 10;
  hl2.FontSize = 10;
  hl3.FontSize = 10;
  fig = gcf;
  fig.Units = 'inches';
  fig.Position(3) = fig.Position(3)*2;
  fig.Position(1) = 1;
  if script_opt.savefig
    figname = fullfile(script_opt.savedir,['LISST_EXPORTSNA_compare_' ptype '_Epoch' num2str(epoch)]);
    standard_printfig_highrespng(figname)
  else
    pause
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