function LISST_plot_EXPORTSNA_EPOCH_compare_sb
if ismac
  base_directory = '/Volumes/toshiba1/LISST_Data_Structured/';
else
  base_directory = 'F:\LISST_Data_Structured\';
end
test_binned = 1;

if test_binned
  DY_sb = fullfile(base_directory,'LISST_sn4025_2021_EXPORTS_DY131','submit','L2','EXPORTS-EXPORTSNA_LISST-Deep_binned_survey_20210504-20210529_R0.sb');
  JC_sb = fullfile(base_directory,'LISST_sn4041_2021_EXPORTS_JC214','submit','L2','EXPORTS-EXPORTSNA_LISST-Deep_binned_process_20210504-20210529_R0.sb');
else
  DY_sb = fullfile(base_directory,'LISST_sn4025_2021_EXPORTS_DY131','submit','L2','EXPORTS-EXPORTSNA_LISST-Deep_survey_20210504-20210529_R0.sb');
  JC_sb = fullfile(base_directory,'LISST_sn4041_2021_EXPORTS_JC214','submit','L2','EXPORTS-EXPORTSNA_LISST-Deep_process_20210504-20210529_R0.sb');
end


[DY, DY_sbHeader, ~] = readsb(DY_sb, 'MakeStructure',true);

[JC, JC_sbHeader, ~] = readsb(JC_sb, 'MakeStructure',true);

keyboard
if ~isequal(DY_sbHeader.fields,JC_sbHeader.fields)
  fprintf('NOT THE SAME FIELDS!\n')
  keyboard
end

fields = strsplit(DY_sbHeader.fields,',');
units  = strsplit(DY_sbHeader.units,',');

DY_table = struct2table(DY);
JC_table = struct2table(JC);

sizes = struct();
sizes.bin_limits = strsplit(DY_sbHeader.psd_bin_size_boundaries,',');
sizes.bins  = nan(32,2);
sizes.width = nan(32,1);
ncnt = 1;
for n = 1:numel(sizes.bin_limits) - 1
  sizes.bins(n,1) = str2double(sizes.bin_limits(ncnt));
  sizes.bins(n,2) = str2double(sizes.bin_limits(ncnt+1));
  sizes.width(n)  = sizes.bins(n,2) - sizes.bins(n,1);
  ncnt = ncnt + 1;
end


%% Get NON differential PSD and VSD
idx_PSD = contains(fields,'PSD_DNSD');
idx_VSD = contains(fields,'PSD_DVSD');
idx_VSF = contains(fields,'VSF');

DYcasts = [19 31 94]; % Disovery found on "Discovery Matchups.xlsx"
JCcasts = [12 23 57]; % JamesCook found on "Cook Matchups.xlsx"
%% Loop through epochs and plot casts
for epoch = 1:3
  %% Plot total volume concentration, mean size, and surface area
  ax = makefig_subplots(3,1); set(gcf,'Name','Exports Epoch Comparison Cast');
  
  for ndata = 1:2
    if ndata == 1
      data   = JC;
      data.cast = str2double(data.station);
      wcast  = find(data.cast == JCcasts(epoch));
      JCwcast = wcast(1);
      clr = 'k';
      CastID = ['JamesCook Cast' num2str(data.cast(wcast(1)))];
      
    elseif ndata == 2
      data   = DY;
      data.cast = str2double(data.station);
      wcast  = find(data.cast == DYcasts(epoch));
      DYwcast = wcast(1);
      clr = 'b';
      CastID = ['Discovery Cast' num2str(data.cast(wcast(1)))];
    end
    
    %% Get NON differential PSD and VSD
    dat_table = struct2table(data);
    
    PSD = table2array(dat_table(wcast,idx_PSD)) .* sizes.width' ;
    VSD = table2array(dat_table(wcast,idx_VSD)) .* sizes.width' ;
    VSF = table2array(dat_table(wcast,idx_VSF));
    
    tot_PSD = sum(PSD,2);
    tot_VSD = sum(VSD,2);
    tot_VSF = sum(VSF,2);
    % Convert from #/m^3 to #/L
    tot_PSD = tot_PSD./1000;
    
    % Convert from uL/m^3 to uL/L
    tot_VSD = tot_VSD./1000;
    
    
    %% Plot profile statistics
    %max_depth = max(data_grid.(depth_var));
    ylims = [0 200]; %max_depth];
    plot(ax(1),data.trans(wcast),data.depth(wcast),'-','Color',clr,'LineWidth',2,'DisplayName',CastID)
    ax(1).XLabel.String = 'Optical transmission [%]';
    hold(ax(1),'on');
    
    plot(ax(2),tot_PSD,data.depth(wcast),'-','Color',clr,'LineWidth',2,'DisplayName',CastID)
    ax(2).XLabel.String = 'Total Particle Concentration [#/L]';
    ax(2).YLabel = []; ax(2).YTickLabel = [];
    hold(ax(2),'on');
    
    plot(ax(3),tot_VSD,data.depth(wcast),'-','Color',clr,'LineWidth',2,'DisplayName',CastID)
    ax(3).XLabel.String = 'Total Volume Concentration [�L/L]';
    ax(3).YLabel = []; ax(3).YTickLabel = [];
    hold(ax(3),'on');
    
    ax(1).Title.String = ['EPOCH ' num2str(epoch)];
    
    ax(1).YLim = ylims;
    ax(2).YLim = ylims;
    ax(3).YLim = ylims;
    ax(1).YDir = 'rev';
    ax(2).YDir = 'rev';
    ax(3).YDir = 'rev';
    
    ax(1).Color = [0.8 0.80 0.8];
    ax(2).Color = [0.8 0.80 0.8];
    ax(3).Color = [0.8 0.80 0.8];
  end
  ax(2).Title.String = {['Discovery Cast '  num2str(DYcasts(epoch)) ' [' num2str(DY.lat(DYwcast),'%.2f') '\circN, ' num2str(DY.lon(DYwcast),'%.2f') '\circW]'];...
                        ['James Cook Cast ' num2str(JCcasts(epoch)) ' [' num2str(JC.lat(JCwcast),'%.2f') '\circN, ' num2str(JC.lon(JCwcast),'%.2f') '\circW]']};
  legend(ax(1),'show','location','sw')
  legend(ax(2),'show','location','se')
  legend(ax(3),'show','location','se')
  
  if test_binned
    figname = fullfile(base_directory,['LISST_Exports_sbfile_binned_compare_Epoch' num2str(epoch)]);
  else
    figname = fullfile(base_directory,['LISST_Exports_sbfile_compare_Epoch' num2str(epoch)]);
  end
  if milliq_compare
    figname = [figname '_milliqcompare'];
  end
  %standard_printfig_highrespng(figname)
  keyboard
end

end