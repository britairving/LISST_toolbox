% dat1 = load('X:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\march10\LISST_sn4025_2018_exports_np_sr1812_binned_L2270456_average300_to_1000m.mat');
% dat2 = load('X:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\march10\LISST_sn4025_2018_exports_np_sr1812_binned_zscat_MilliQ_avg.mat');
% dat3 = load('X:\LISST_Data\LISST_sn4025_2018_exports_np_sr1812\march10\LISST_sn4025_2018_exports_np_sr1812_binned_zscat_FSW_avg.mat');

%% 0 | Set up script options
opt = dat1.opt; % all the same - just use first
depth_var = 'depth';% 'depth_orig'; % 'depth'
date_var  = 'date'; % 'date_orig';  % 'date'
script_opt.axcolor = [0.8 0.80 0.8];
script_opt.savefig = 0;
casts = dat1.data_grid.cast;
% casts = [10, 68];

num_casts = numel(casts);
max_depth = max(dat1.data_grid.(depth_var));
dias      = dat1.dat.file1.dias;     % dat.meta.dias = 'midpoint of size bins'
dep       = dat1.data_grid.(depth_var);

data_grid1 = dat1.data_grid; % InSitu
data_grid2 = dat2.data_grid; % MilliQ
data_grid3 = dat3.data_grid; % FSW
[~,zname1,~] = fileparts(dat1.cfg.proc_options.zscfile{1});
[~,zname2,~] = fileparts(dat2.cfg.proc_options.zscfile{1});
[~,zname3,~] = fileparts(dat3.cfg.proc_options.zscfile{1});

%% scattering 
% scat(i,:)=data(i,1:32)/tau(i)-(zsc(1:32)*data(i,36)/zsc(36));

for nc = 1:num_casts
  wcast = dat1.data_grid.cast == casts(nc);
  CastID = ['Cast' num2str(dat1.data_grid.cast(wcast))];
  % Load BRU & DAT files
  fprintf(' Plotting LISST %s data\n',CastID)

  %% Plot total volume concentration, mean size, and surface area
  ax = makefig_subplots(3,1); set(gcf,'Name',dat1.opt.project);
  
  for ndata = 1:3
    if ndata == 1
      data_grid = data_grid1;
      clr = 'b';
      lstring = zname1;
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
    elseif ndata == 2
      data_grid = data_grid2;
      clr = 'k';
      lstring = zname2;
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
    elseif ndata == 3
      data_grid = data_grid3;
      clr = 'r';
      lstring = zname3;
      [~,datfile,~] = fileparts(data_grid.datfile{wcast});
      [~,zscfile,~] = fileparts(data_grid.zscfile{wcast});
    end
    lstring = strrep(lstring,'_','\_');
    %% Plot profile statistics
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
    ax(2).Title.String = {[strrep(CastID,'_', ' ') ': ' datfile];['zscat: ' strrep(zscfile,'_', '\_')]};
    hold(ax(2),'on');
    
    x = data_grid.surf_area(wcast,:); % squeeze(data_grid.cscat(wcast,:,:));
    igood = isfinite(x); x = x(igood); y = dep(igood);
    plot_against_pressure(ax(3),x,y,ylims,clr,lstring);
    ax(3).XLabel.String = 'Surface Area [cm^2/l]';
    ax(3).YLabel = []; ax(3).YTickLabel = [];  ax(3).Title.String = [];
    hold(ax(3),'on');
  end
  legend(ax(1),'show','location','se')
  legend(ax(2),'show','location','se')
  legend(ax(3),'show','location','se')
  if script_opt.savefig
    figname = fullfile(opt.path.figs,[opt.project '_' CastID '_zscat_compare']);
    standard_printfig_lowrespng(figname)
  else
    keyboard
  end
end
keyboard
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
