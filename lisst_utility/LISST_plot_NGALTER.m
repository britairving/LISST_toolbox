function LISST_plot_NGALTER(l2_fname)
%FUNCTION LISST_plot_NGALTER
%
%  Syntax:
%     LISST_plot_NGALTER(l2_fname)
%  
%  Description:
%     Simple function that reads in LISST data formatted for NGA LTER
%     level2 submission to AXIOM. 
%
%  Inputs: 
%     l2_fname | filename with path to L2 csv file 
%
%  Authors:
%     Brita K Irving  <bkirving@alaska.edu>
%% 
close all;
ignore_flagged = 1; % 1 = removes flagged data before plotting, 0 = plots everything, regardless of flags
save_figure    = 1; % 1 = saves figures, 0 = does not.
[flder,name,~] = fileparts(l2_fname);

%% Read column names and units
fileID = fopen(l2_fname,'r');
fline1 = fgetl(fileID);
fline2 = fgetl(fileID);
fclose(fileID);
% convert to cell array
names = strsplit(fline1,',');
units = strsplit(fline2,',');
% parse bin size_bin from column name
size_bin = struct();
size_bin.sizes = [];
size_bin.names = {};
size_bin.med   = [];

if any(contains(names,'PSD'))
  sPSD = names(contains(names,'PSD') & contains(names,'#/L')); % only pull out PSD [#/L], not PSD_DNSD or others
  for ns = 1:numel(sPSD)
    str = strsplit(sPSD{ns},'_');
    szs = strsplit(str{2},'-');
    size_bin.names{ns} = str{2};
    size_bin.sizes(ns,1) = str2double(szs{1});
    size_bin.sizes(ns,2) = str2double(erase(szs{2},'um'));
    % Calculate middle of bin as the 'ssquare root of the product of the lower and upper bin limit
    size_bin.med(ns) = sqrt(size_bin.sizes(ns,1)*size_bin.sizes(ns,2));
  end
end 

%% Read data
topt = detectImportOptions(l2_fname);
topt.VariableUnitsLine = 2;
topt.DataLines(1) = 3;
data = readtable(l2_fname,topt);
% Examine data
% head(data);

%% Convert to table so easier to look at and plot
idx_PSD = contains(names,'PSD') & contains(names,'#/L');
idx_VSD = contains(names,'VSD') & contains(names,'uL/L');
PSD = table2array(data(:,idx_PSD));
VSD = table2array(data(:,idx_VSD));
data(:,idx_PSD | idx_VSD) = [];
data = table2struct(data,'ToScalar',true);
data.PSD = PSD;
data.VSD = VSD;

%% Loop through casts and plot data vs depth
% Plots data vs depth for variables listed in vars, generates a single plot
% for entire cruise.
vars = {'transmission____' 'mean_size__um_' 'PSD' 'VSD'};
ax = makefig_subplots(numel(vars),1);
casts = unique(data.Cast_Number);
for nc = 1:numel(casts)
  if ignore_flagged
    idx_cast = data.Cast_Number == casts(nc) & data.quality_flag <= 2; % not_evaulated or good
  else
    idx_cast = data.Cast_Number == casts(nc);
  end
  
  for na = 1:numel(vars)
    % average and smooth the data
    dat = nanmean(data.(vars{na})(idx_cast,:),2); % if this is an Nx1 variable, taking the mean will do nothing, if this is a NxM variable (e.g. PSD or VSD), it will average all size bins
    dat = smoothdata(dat,'movmedian',125); % not binned, so smooth a lot....
    
    plot(ax(na),dat,data.Depth__m_(idx_cast),'-');%,'DisplayName',['Cast' num2str(casts(nc))])
    
    ax(na).YDir = 'reverse';
    ax(na).YLim(1) = 0;
    ax(na).LineWidth = 3;
    if na ~=1
      ax(na).YTickLabel = [];
    end
    if nc == 1
      ax(4).XScale = 'log'; % VSD
      ax(3).XScale = 'log'; % PSD
      %ax(4).XScale = 'linear'; % VSD
      %ax(3).XScale = 'linear'; % PSD
      hold(ax(na),'on');
      %grid(ax(na),'on');
    end
  end % LOOP THROUGH AXES
end % LOOP THROUGH CASTS
ax(1).YLabel.String = 'Depth [m]';
ax(1).XLabel.String = 'transmission [%]';
ax(2).XLabel.String = 'Mean Size [\mum]';
ax(3).XLabel.String = 'Average PSD [#/L]';
ax(4).XLabel.String = 'Average VSD [\mul/l]';
ax(2).Title.String  = [strrep(name,'_','\_') ' all profiles vs depth'];
if save_figure
  try
    standard_printfig_highrespng(fullfile(flder,name));
  catch
    print(fullfile(flder,name),'-dpng','-r250');
  end
end

%% Plot PSD and VSD for each size class
cmap = lansey(32);
ax = makefig_subplots(2,1);
set(gcf,'Name',[name ' VSD and PSD by depth'])
casts = unique(data.Cast_Number);
for nc = 1:numel(casts)
  if ignore_flagged
    idx_cast = find( data.Cast_Number == casts(nc) & data.quality_flag <= 2); % not_evaulated or good
  else
    idx_cast = find( data.Cast_Number == casts(nc));
  end
  scast = ['Cast' num2str(casts(nc))];
  if isempty(idx_cast)
    continue
  else
    %% PSD and VSD for each size class at all depths within one station
    for n = 1:32
      try
        cPSD = smoothdata(squeeze(data.PSD(idx_cast,n)),'movmedian',125);
        cVSD = smoothdata(squeeze(data.VSD(idx_cast,n)),'movmedian',125);
        
        hPSD(n) = plot(ax(1), cPSD,data.Depth__m_(idx_cast),'Color',cmap(n,:),'Marker','s','LineStyle','-','LineWidth',2,'DisplayName',size_bin.names{n});
        hVSD(n) = plot(ax(2), cVSD,data.Depth__m_(idx_cast),'Color',cmap(n,:),'Marker','s','LineStyle','-','LineWidth',2,'DisplayName',size_bin.names{n});
      catch
        keyboard
      end
      if n == 1
        hold(ax(1),'on');hold(ax(2),'on');
        %grid(ax(1),'on');grid(ax(2),'on');
      end
    end
    %ax(1).XScale = 'linear';
    %ax(2).XScale = 'linear';
    ax(1).XScale = 'log';
    ax(2).XScale = 'log';
    ax(1).YDir   = 'rev';
    ax(2).YDir   = 'rev';
    
    ax(1).XLabel.String = 'PSD [#/L]';
    ax(2).XLabel.String = 'VSD [\muL/L]';
    ax(1).YLim(1) = 0;
    ax(2).YLim(1) = 0;
    hl = legend(ax(2),'show'); hl.FontSize = 14;
    hl.Location = 'bestoutside';
    ylabel(ax(1),'Depth [m]');
    try
      title(ax(1),[data.Station{idx_cast(1)} ' | ' scast]);
    catch
      keyboard
    end
    % Save figure
    if save_figure
      figname = fullfile(flder,[scast '_' data.Station{idx_cast(1)} '.png']);
      try
        standard_printfig_highrespng(figname);
      catch
        print(figname,'-dpng','-r250');
      end
    else
      fprintf('paused.. press any key to continue\n')
      pause
    end
        
    % Clear axes
    cla(ax(1));
    cla(ax(2));
  end
  
end

end %% MAIN FUNCTION 