function [cfg, downcasts] = LISST_select_background_scatterfiles(cfg,downcasts)
%function LISST_select_background_scatterfiles
%
%  Syntax:
%    [cfg,downcasts] = LISST_select_background_scatterfiles(cfg,downcasts)
%
%  Description: 
%     
%  Justification: 
%    The background (also called zscat) MUST BE measured prior to each
%    experiment. The field data is the summation of the background and
%    contributions from particles. The background as well as the scattering
%    from particles, both is attenuated in water. Consequently, the data
%    must be de-attenuated and for this purpose, the transmission
%    measurement is included in the data stream.
%
%  REQUIRED INPUT: 
%    Section 2 requires user input to select background scatterfiles for
%    each .dat file.
% 
%    Can use script LISST_select_background_scatterfiles.m to visualize
%    available background files ... 
%  
%    IF no specific information is added for the project, an average
%    background scatterfile is created using the average of all ascii files
%    in the [cfg.project]\background_scatter_files\ folder. 
%
%  Refereces:
%    http://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/
%    LISST-Deep-Users-Manual-May-2013.pdf
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Set script defaults
% Assume a max of 4 water types 
% FSW = filtered seawater
% MQW = MilliQ water
% DIQ = De-ionized water
% UNK = unknown
cmaps  = {'cool' 'copper' 'bone' 'autumn'}; % colormaps used to see differences in backgrounds
symbls = {'p-'     'o-'     '<--'   'x-'};     % symbols used to see differences in backgrounds
fprintf(' \n------------------------------------\n')
fprintf(' LISST_select_background_scatterfiles\n')
%% 1 | Load background scatter files
zscat = struct();
% Load factory background
[~,factory_zscat,ext] = fileparts(cfg.path.file_factoryz);
zscat.factory.file = [factory_zscat ext];
zscat.factory.data = load(cfg.path.file_factoryz);
zscat.factory.date = datenumfromdata(cfg.year,zscat.factory.data(39:40)');
% store factory laser reference
cfg.inst.flref = zscat.factory.data(36); %
% Store all backgrounds to single matrix so can pull out max and min for
% uncertainty estimation
znm_all = {zscat.factory.file};
zsc_all = zscat.factory.data;


if isfield(cfg,'background')
  % initialize variable
  rm_water_type = {};
  water_types = fieldnames(cfg.background);
  for ntype = 1:numel(water_types)
    wt = water_types{ntype};
    
    if isempty(cfg.background.(wt))
      % These can be empty because all types are kept in the
      % LISST_project_config.m script to serve as a template
      rm_water_type = [rm_water_type; wt];
    else
      zscat.(wt) = struct();
      zscat.(wt).file = cfg.background.(wt);
      % initialize background measurement array to store all background
      % data plus mean & median
      zscat.(wt).data = nan(40,numel(zscat.(wt).file));
      % Loop through background scatter files and load the data
      for nz = 1:numel(cfg.background.(wt))
        zfile = fullfile(cfg.path.dir_zsc,zscat.(wt).file{nz});
        try
          zscat.(wt).data(:,nz) = load(zfile);
          znm_all = [znm_all, zscat.(wt).file(nz)];
          zsc_all = [zsc_all, zscat.(wt).data(:,nz)];
        catch
          fprintf(' could not load file %s\n',zfile);
        end
      end
      zscat.(wt).date = datenumfromdata(cfg.year,zscat.(wt).data(39:40,:)')'; % 39th element = (Day*100 + Hour) at which data taken
      % remove copies and sort by date 
      [~,isort] = unique(zscat.(wt).date); 
      zscat.(wt).data = zscat.(wt).data(:,isort);
      zscat.(wt).file = zscat.(wt).file(isort);
      zscat.(wt).date = zscat.(wt).date(isort);
      % Calculate the median and the average
      avg_data = mean(zscat.(wt).data,2);   % mean background for water type
      avg_date = mean(zscat.(wt).date);
      med_data = median(zscat.(wt).data,2); % median background for water type
      med_date = median(zscat.(wt).date);
      % add to mean and median 
      zscat.(wt).file = [zscat.(wt).file, [wt ' average'], [wt ' median']];
      zscat.(wt).date = [zscat.(wt).date, avg_date, med_date];
      zscat.(wt).data = [zscat.(wt).data, avg_data, med_data];
    end
   
  end
else
  fprintf(' Need to set this up!\n')
  keyboard
end
% Remove unnecessary water_types
water_types(contains(water_types,rm_water_type)) = [];

%% Pull out minimum and maximum zscat for uncertainty calculation
if cfg.proc_options.calculate_range
  cfg.range.header   = 'Range is the difference of the parameter calculated with the maximum background and the parameter calculated with the minimum background';
  cfg.range.range1.zscfile = {}; % range1 = difference of data (lowest zscat that is NOT the factory) and data (highest zscat)
  cfg.range.range2.zscfile = {}; % range2 = difference of data (factory zscat) and data (highest zscat). This will be larger than range1.
  
  % remove factory backgrounds
  zsc_all(:,1) = [];
  znm_all(1)   = [];
  
  % Pull out minimum and maximum backgrounds - taken as the sum of 
  [~,idx_minzscat] = min(sum(zsc_all(1:32,:)));
  [~,idx_maxzscat] = max(sum(zsc_all(1:32,:)));
  
  cfg.range.range1.zscfile = {znm_all{idx_minzscat} znm_all{idx_maxzscat}};
  cfg.range.range2.zscfile = {zscat.factory.file    znm_all{idx_maxzscat}};
  clearvars znm_all zsc_all
  
end

%% Calculate transmission with different backgrounds
% Transmission t is computed from the ratio of laser power transmitted
% through water in turbid conditions to its value in clean conditions –
% this is the ratio of variable 33 measured in situ and its value in clean
% water. Since laser output can drift over time, a correction for drift is
% applied using the laser reference sensor. If r is the ratio of laser
% transmitted power to laser reference measurement in clean water, then it
% follows that despite any drift in laser output,  the clean water
% transmitted power would be the product of laser reference and r, i.e. r
% *element(36). It then follows that transmission for the ith data record
% in-situ is
% First, calculate the transmission with the factory background
r = zscat.factory.data(33)/zscat.factory.data(36);
zscat.factory.tau = LISST_calculate_transmission(cfg,downcasts,r);

for ntype = 1:numel(water_types)
  wt = water_types{ntype};
  fprintf('  Calculating transmission using %s backgrounds\n',wt)
  % compute the laser power/laser reference ratio (to adjust for drift in laser output power over time)
  for nz = 1:numel(zscat.(wt).date)
    r = zscat.(wt).data(33,nz)/zscat.(wt).data(36,nz);
    tau = LISST_calculate_transmission(cfg,downcasts,r);
    zscat.(wt).tau(nz) = tau;
  end % Loop through background files for water type
end % Loop through water types


%% Calculate in-situ background
zscat = LISST_select_insitu_background(cfg,downcasts,zscat);

%% 2 | Plot all background data
fprintf(' Plotting available background scatterfiles\n')
% Initialize plot handles
h1 = struct();
h2 = struct();
h3 = struct();

% Create figure and plot factory
makefig; 
makefig; ax1 = gca;  % Backscatter value per bin
makefig; ax2 = gca;  % tau vs time
makefig; ax3 = gca;  % tau vs depth

% Plot factor backscatter
h1.factory = plot(ax1,1:32,zscat.factory.data(1:32),'k-','LineWidth',3,'DisplayName',['Factory ' datestr(zscat.factory.date,'yyyy-mm-dd')]);
h2.factory = plot(ax2,zscat.factory.tau.time_m,zscat.factory.tau.surf,'k-','LineWidth',3,'DisplayName',['Factory ' datestr(zscat.factory.date,'yyyy-mm-dd')]);
h3.factory = plot(ax3,zscat.factory.tau.mean,zscat.factory.tau.depth_m,'k-','LineWidth',3,'DisplayName',['Factory ' datestr(zscat.factory.date,'yyyy-mm-dd')]);

hold(ax1,'on'); grid(ax1,'on');
hold(ax2,'on'); grid(ax2,'on');
hold(ax3,'on'); grid(ax3,'on');

if isfield(zscat,'insitu')
  % Plot factor backscatter
  h1.insitu = plot(ax1,1:32,zscat.insitu.data(1:32),'b-','LineWidth',3,'DisplayName',['insitu ' zscat.insitu.file]);
  h2.insitu = plot(ax2,zscat.insitu.tau.time_m,zscat.insitu.tau.surf,'b-','LineWidth',3,'DisplayName',['in-situ ' zscat.insitu.file]);
  h3.insitu = plot(ax3,zscat.insitu.tau.mean,zscat.insitu.tau.depth_m,'b-','LineWidth',3,'DisplayName',['in-situ ' zscat.insitu.file]);
end

% Title
ax1.Title.String = strrep(cfg.project,'_','\_');
ax2.Title.String = strrep(cfg.project,'_','\_');
ax3.Title.String = strrep(cfg.project,'_','\_');


% transmission vs depth
ax3.XLabel.String = 'Transmission';
ax3.YLabel.String = 'Depth [m]';
ax3.YDir = 'reverse';

% transmission vs time
datetick(ax2,'x','mm/dd','keeplimits','keepticks')
ax2.YLabel.String = 'Transmission averaged from 0-200m';

% backscatter vs bin
ax1.XLabel.String = 'Bin'; 
ax1.YLabel.String = 'Backscatter Value';

% plot background measurements 
nzscat = 0; % Total number of backgrounds 
for nw = 1:numel(water_types)
  % Pull out water_type name
  wt   = water_types{nw}; 
  clrs = eval([cmaps{nw} '(' num2str(numel(zscat.(wt).date)) ')']);
  
  % Loop through all backgrounds taken with this water type
  for nz = 1:numel(zscat.(wt).date)
    nzscat = nzscat + 1;
    if contains(zscat.(wt).file(nz),{'average' 'median'})
      lnwdth = 2.0;
    else
      lnwdth = 1.0;
    end
    % Backscatter value per bin
    h1.(wt)(nz) = plot(ax1, 1:32,zscat.(wt).data(1:32,nz),symbls{nw},...
                       'Color',clrs(nz,:),'LineWidth',lnwdth,'MarkerFaceColor',clrs(nz,:),'MarkerSize',8,...
                       'DisplayName',['<' num2str(nzscat) '> ' wt ' ' zscat.(wt).file{nz}]);
    % tau vs time
    h2.(wt)(nz) = plot(ax2,zscat.(wt).tau(nz).time_m,zscat.(wt).tau(nz).surf,symbls{nw},...
                       'Color',clrs(nz,:),'LineWidth',lnwdth,'MarkerFaceColor',clrs(nz,:),'MarkerSize',4,...
                       'DisplayName',['<' num2str(nzscat) '> ' wt ' ' zscat.(wt).file{nz}]); 
    % tau vs depth
    h3.(wt)(nz) = plot(ax3,zscat.(wt).tau(nz).mean,zscat.(wt).tau(nz).depth_m,symbls{nw},...
                       'Color',clrs(nz,:),'LineWidth',lnwdth,'MarkerFaceColor',clrs(nz,:),'MarkerSize',4,...
                       'DisplayName',['<' num2str(nzscat) '> ' wt ' ' zscat.(wt).file{nz}]);    
  end % Loop through all backgrounds taken with this water type
end % Loop through all water types

%% legends
for a = 1:3
  switch a
    case 1 % backscatter value per bin
      st = 'backscatter_per_bin';
      ax = ax1;
    case 2 % tau vs time
      st = 'tau_vs_time';
      ax = ax2; 
      plot_transmission_flags(ax,'y'); % plot the flags
    case 3 % tau vs depth
      st = 'tau_vs_depth';
      ax = ax3; 
      plot_transmission_flags(ax,'x'); % plot the flags
  end
  
  try
    hl = legend(ax,'show');
    hl.FontSize = 14;
    hl.Location = 'northeastoutside'; %hl.Position = [0.75 0.42 0.23 0.57]; % move to upper right corner
    hl.Interpreter = 'none' ; % don't have to use strrep(text,'_','\_')
  catch % R2019a running on Mac has issues with legend..
    %ht = legend_alternative(ax1);
    ht = text(ax);
    try
      hl = legend(ax,'show');
      hl.FontSize = 14;
      hl.Location = 'northeastoutside'; %hl.Position = [0.75 0.42 0.23 0.57]; % move to upper right corner
      hl.Interpreter = 'none' ; % don't have to use strrep(text,'_','\_')
      delete(ht);
    catch
      ht = legend_alternative(ax);
    end
    
  end
  %% Save Figure
  if cfg.savefig
    figure(ax.Parent);
    axes(ax);
    figname = fullfile(cfg.path.dir_figs,[cfg.project '_' st]);
    standard_printfig_highrespng(figname)
  end
end

if strcmp(datestr(now,'yyyymmdd'),'20210324')
  % TEST MULTIPLE BACKGROUNDS
  % save median FSW
  zscat.FSW.file{end} = 'FSW_zsc_median.asc';
  zscat_file = fullfile(cfg.path.dir_zsc,zscat.FSW.file{end});
  fprintf(' Saving insitu background to %s\n', zscat_file)
  fileID = fopen(zscat_file,'w');
  fprintf(fileID,'%0.4f\n',zscat.FSW.data(:,end));
  fclose(fileID);
  % save median MilliQ
  zscat.MQW.file{end} = 'MilliQ_zsc_median.asc';
  zscat_file = fullfile(cfg.path.dir_zsc,zscat.MQW.file{end});
  fprintf(' Saving insitu background to %s\n', zscat_file)
  fileID = fopen(zscat_file,'w');
  fprintf(fileID,'%0.4f\n',zscat.MQW.data(:,end));
  fclose(fileID);
  
  % PROCESS 4 ways to compare...
  downcasts.zscfile = [cellstr(repmat(zscat.factory.file,size(downcasts.datfile))),... % 1. factory background
    cellstr(repmat(zscat.insitu.file,size(downcasts.datfile))),...  % 2. insitu
    cellstr(repmat(zscat.FSW.file{end},size(downcasts.datfile))),...% 3. median FSW
    cellstr(repmat(zscat.MQW.file{end},size(downcasts.datfile)))];  % 4. median MilliQ
  
  cfg.proc_options.zscat_choice  = {'factory_zscat', 'insitu_zscat', 'median_FSW_zscat', 'median_MilliQ_zscat'};
  return
end

%% Select which backgrounds to use...
% Bring background vs bin plot to the front
axes(ax1)
zscat_types = fieldnames(zscat);
fprintf(' \n')
fprintf(' ---------------------------------------------------\n')
fprintf(' Which type of background do you want to use for processing?\n');
done = 0; 
while ~done
  for nw = 1:numel(zscat_types)
    fprintf('   <%d>  %s\n',nw,zscat_types{nw})
  end
  fprintf('   <99>  STOP\n')
  chc_type = input('  Enter choice: ');
  if chc_type == 99
    fprintf(' stopped here\n')
    keyboard
  else
    wt = zscat_types{chc_type};
    zscat = zscat.(wt);
    done  = 1;
  end
end %% WHILE ~DONE
% Delete other plots
rm_these = find(~contains(zscat_types,zscat_types{chc_type}));
for nrm = 1:numel(rm_these)
  rmtype = zscat_types{rm_these(nrm)};
  delete(h1.(rmtype));
  delete(h2.(rmtype));
  delete(h3.(rmtype));
  h1 = rmfield(h1,rmtype);
  h2 = rmfield(h2,rmtype);
  h3 = rmfield(h3,rmtype);
end
% Renumber legend
try pause(1);end
remove_plots = [];
for np = 1:numel(h1.(wt))
  if contains(h1.(wt)(np).DisplayName,{'median' 'average'})
    remove_plots = [remove_plots; np];
  end
  idx  = strfind(h1.(wt)(np).DisplayName,'>');
  lstr = h1.(wt)(np).DisplayName(1:idx);
  h1.(wt)(np).DisplayName = strrep(h1.(wt)(np).DisplayName,lstr,['<' num2str(np) '>']);
  h2.(wt)(np).DisplayName = strrep(h2.(wt)(np).DisplayName,lstr,['<' num2str(np) '>']);
  h3.(wt)(np).DisplayName = strrep(h3.(wt)(np).DisplayName,lstr,['<' num2str(np) '>']);
end

delete(h1.(wt)(remove_plots));
delete(h2.(wt)(remove_plots));
delete(h3.(wt)(remove_plots));
% reset axis limits
try
  ax2.YLim(1) = min(nanmin([zscat.tau.surf]));
  ax3.XLim(1) = min(nanmin([zscat.tau.surf]));
end
%% Visually highlight selected background measurements
if numel(zscat.file) > 1
  
  fprintf(' Do you want to highlight any of the background measurements on the plot?\n')
  chc1 = input('  Enter choice <1/0> ');
  if chc1
    while chc1
      for np = 1:numel(h1.(wt))
        if ~isvalid(h1.(wt)(np))
          continue
        end
        fprintf('   %s\n',h1.(wt)(np).DisplayName)
      end
      chc = input('  Enter choice (ret to skip): ');
      if ~isempty(chc)
        orig_lw  = h1.(wt)(chc).LineWidth;
        orig_clr = h1.(wt)(chc).Color;
        h1.(wt)(chc).LineWidth = 5;
        h1.(wt)(chc).Color = 'r';
        uistack(h1.(wt)(chc),'top');
        chc1 = input('Do you want to highlight anymore of the background measurements? <1/0> ');
        % reset to original width and color
        h1.(wt)(chc).LineWidth = orig_lw;
        h1.(wt)(chc).Color     = orig_clr;
      else
        chc1 = 0;
      end
    end
    chc1 = 1;
  end % highlight selected backscatter files on plot
end


%% Remove files
fprintf(' \n')
fprintf(' Do you want to remove any backgrounds from consideration?\n');
chc_rm = input('  Enter choice <1/0>: ');
if chc_rm
  remove_zscat = 1;
  remove_these = []; % store indices of background files selected to remove
  while remove_zscat
    for nz = 1:numel(h1.(wt))
      % If handle was delted it will still be a variable, but will be
      % invalid
      if ~isvalid(h1.(wt)(nz))
        continue
      end
      if ismember(nz,remove_these)
       fprintf('   <%d>\t----------------\n',nz)
      else
        fprintf('   %s\n',h1.(wt)(nz).DisplayName)
      end
    end
    idx_rm = input('  Enter choice: ');
    ht = ['htype' num2str(chc_type)];
    % highlight select file on plot    
    orig_lw  = h1.(wt)(idx_rm).LineWidth;
    orig_clr = h1.(wt)(idx_rm).Color;
    h1.(wt)(idx_rm).LineWidth = 5;
    h1.(wt)(idx_rm).Color     = 'r';
    fprintf(' \n')
    fprintf('   You selected to ignore %s\n',h1.(wt)(idx_rm).DisplayName)
    chc_rm = input('  Is that correct? <1/0> ');

    if chc_rm
      remove_these = [remove_these; idx_rm];
      h1.(wt)(idx_rm).Visible = 'off';
    else
      % reset to original linewidth and color
      h1.(wt)(idx_rm).LineWidth = orig_lw;
      h1.(wt)(idx_rm).Color     = orig_clr;
    end
    fprintf(' \n')
    chc_again = input('  Do you want to remove any other files? <1/0> ');
    if isempty(chc_again) || chc_again == 0
      remove_zscat = 0; % exit from loop
    end
  end % WHILE remove_files = loop through and remove select files until finished
  % Remove background(s) from all zscat fields
  zscat.file(remove_these)   = [];
  zscat.date(remove_these)   = [];
  zscat.data(:,remove_these) = [];
  zscat.tau(remove_these)    = [];
  % Remove plot(s) from all figures
  delete(h1.(wt)(remove_these));
  delete(h2.(wt)(remove_these));
  delete(h3.(wt)(remove_these));
  
  % save figure
  if cfg.savefig
    figname = fullfile(cfg.path.dir_figs,[cfg.project '_background_scatter_measurements_trimmed']);
    standard_printfig_lowrespng(figname)
  end
end % opted to remove some of the background files (not deleting them, just removing from plot and average/interpolating/nearest below)

%% Select which method to assign zscat files to process dat files (mean, median, individual, or interpolated)
% pull out dates for each cast
[casts,iu] = unique(downcasts.cast,'stable');
dates= downcasts.date(iu);
% Remove mean and median background in case any were removed
rm_idx = contains(zscat.file,{' average' ' median'}); % these should be the last two, but do this just in case
zscat.file(rm_idx)   = [];
zscat.data(:,rm_idx) = [];
zscat.date(rm_idx)   = [];
zscat = rmfield(zscat,'tau'); % don't need this anymore, so avoid confusion by removing it.
% reset downcasts.zscfile
downcasts.zscfile = repmat({''},size(downcasts.datfile));
done_method = 0;
while ~done_method
  fprintf(' \n')
  fprintf(' ---------------------------------------------------\n')
  fprintf(' Type of backgrounds selected: %s\n',wt);
  fprintf(' ---------------------------------------------------\n')
  fprintf(' Select which method to apply backgrounds for each datfile\n');
  fprintf('   <1>   Nearest background measurements in time\n')               % Assign background to each dat file using nearest time
  fprintf('   <2>   Interpolated background measurements\n')                  % Interpolate background measurements using time
  fprintf('   <3>   Mean of background measurements across whole cruise\n')   % Mean already calculated, but recalculate if any zscat files were removed
  fprintf('   <4>   Median of background measurements across whole cruise\n') % Median already calculated, but recalculate if any zscat files were removed
  fprintf('   <5>   Single background file\n')                                % Use a single background file to process all data
  fprintf('   <99>  Stop\n')
  chc_method = input('  Enter choice: ');
  switch chc_method
    %% Use background measurements closest in time to start of datfile
    case 1 % NEAREST IN TIME
      done_method = 1;
      cfg.proc_options.zscat_choice = ['nearest' wt 'zscat'];
      for  nf = 1:numel(casts)
        % Pull out downcasts index for this cast
        idx = downcasts.cast == casts(nf);
        % find index of closest background measurement in time
        nr = nearest(zscat.date,dates(nf));
        cfg.proc_options.zscfile{nf}        = zscat.file{nr};
        downcasts.zscfile(idx) = cfg.proc_options.zscfile(nf);
      end
    %% Use interpolated background measurements for each dat file
    case 2 % INTERPOLATE
      cfg.proc_options.zscat_choice = ['interp' wt 'zscat'];
      % Initialize interpolated zscat structure
      zsc_int  = nan(40,numel(casts));
      for nl = 1:40
        zsc_int(nl,:) = interp1(zscat.date,zscat.data(nl,:),dates);
        % interpolate leading and trailing values too
        if any(any(isnan(zsc_int(nl,:))))
          zsc_int(nl,:) = fillmissing(zsc_int(nl,:)','linear','EndValues','extrap')';
        end
      end

      % create directory
      cfg.path.dir_zsc = fullfile(cfg.path.dir_zsc,['interpolated_' wt '_' datestr(now,'yyyymmdd_HHMM')]);
      if ~exist(cfg.path.dir_zsc,'dir')
        fprintf('  making directory %s\n',cfg.path.dir_zsc)
        mkdir(cfg.path.dir_zsc)
      end
      % save interpolated backgrounds to files
      for nf = 1:numel(casts)
        % Pull out downcasts index for this cast
        idx = find(downcasts.cast == casts(nf));
        [~,datname,~] = fileparts(downcasts.datfile{idx(1)});
        zscatname = [datname '_interpolated_zscat_' wt '.asc'];
        fullzscatname = fullfile(cfg.path.dir_zsc,zscatname);
        fprintf(' writing %s\n',fullzscatname)
        fileID = fopen(fullzscatname,'w');
        fprintf(fileID,'%0.4f\n',zsc_int(:,nf));
        fclose(fileID);
        %save to file cfg structure
        cfg.proc_options.zscfile{nf}        = zscatname;
        downcasts.zscfile(idx) = cfg.proc_options.zscfile(nf);
      end
      done_method = 1;
    %% Use averaged background measurements over whole cruise
    case 3 % MEAN
      cfg.proc_options.zscat_choice = ['mean' wt 'zscat'];
      avg_zscat = mean(zscat.data,2); % average background for water type
      
      zscatname = [cfg.project '_mean_' wt '_zscat.asc'];
      cfg.proc_options.zscfile(:)       = {zscatname};
      downcasts.zscfile(:) = {zscatname};
      
      fprintf(' \nAverage background calculated from files: %s\n',strjoin(zscat.file,', '))
      fprintf(' \nSaving averaged background scatter measurements to file: %s\n',zscatname)
      fileID = fopen(fullfile(cfg.path.dir_zsc,zscatname),'w');
      fprintf(fileID,'%.3f\n',avg_zscat);
      fclose(fileID);
      done_method = 1;
      % update figure
      plot(ax1,1:32,avg_zscat(1:32),'k-','LineWidth',5,'DisplayName',['SELECTED: mean zscat from ' wt]);
      % save figure
      if cfg.savefig
        figname = fullfile(cfg.path.dir_figs,[cfg.project '_background_scatter_measurements_selected']);
        standard_printfig_lowrespng(figname)
      end
    %% Use median of background measurements over whole cruise
    case 4 % MEDIAN
      cfg.proc_options.zscat_choice = ['median' wt 'zscat'];
      med_zscat = median(zscat.data,2); % average background for water type
      zscatname = [cfg.project '_median_' wt '_zscat.asc'];
      cfg.proc_options.zscfile(:)       = {zscatname};
      downcasts.zscfile(:) = {zscatname};
      fprintf(' \nMedian background calculated from files: %s\n',strjoin(zscat.file,', '))
      fprintf(' \nSaving median background scatter measurements to file: %s\n',zscatname)
      fileID = fopen(fullfile(cfg.path.dir_zsc,zscatname),'w');
      fprintf(fileID,'%.3f\n',med_zscat);
      fclose(fileID);
      done_method = 1;
      % update figure
      plot(ax1,1:32,med_zscat(1:32),'k-','LineWidth',5,'DisplayName',['SELECTED: median zscat from ' wt]);
      % save figure
      if cfg.savefig
        figname = fullfile(cfg.path.dir_figs,[cfg.project '_background_scatter_measurements_selected']);
        standard_printfig_lowrespng(figname)
      end
    %% Use single background file for all data
    case 5 % SINGLE FILE
      cfg.proc_options.zscat_choice = ['single' wt 'zscat'];
      done_single_file = 0;
      if numel(zscat.file) == 1
        zscatname = zscat.file{1};
        cfg.proc_options.zscfile(:)       = {zscatname};
        downcasts.zscfile(:) = {zscatname};
        done_single_file = 1;
      else
        while ~done_single_file
          for nz = 1:numel(zscat.file)
            fprintf('  <%d> %s\n',nz,zscat.file{nz})
          end
          chc_file = input(' Enter choice: ');
          try
            zscatname = zscat.file{chc_file};
            cfg.proc_options.zscfile(:)       = {zscatname};
            downcasts.zscfile(:) = {zscatname};
            done_single_file = 1;
          catch
            fprintf(' could not assign choice\n')
            done_single_file = 0;
          end
        end
        % update figure
        h1.(wt)(chc_file).LineWidth = 5;
        h1.(wt)(chc_file).Color = 'k';
        h1.(wt)(chc_file).DisplayName = ['Selected: ' h1.(wt)(chc_file).DisplayName];
      end

      % save figure
      if cfg.savefig
        figname = fullfile(cfg.path.dir_figs,[cfg.project '_background_scatter_measurements_selected']);
        standard_printfig_lowrespng(figname)
      end
      done_method = 1;
    %% STOP HERE
    case 99
      fprintf(' Entering keyboard mode.. enter "dbcont" to continue\n')
      keyboard
    otherwise
      fprintf(' Incorrect selection... try again\n')
      done_method = 0;
  end
end %% Select which method to use to assign background file for processing datfile

%% FUNCTION PLOT_TRANSMISSION_FLAGS
  function plot_transmission_flags(ax,which_ax)
    ylim = ax.YLim; % store original yaxis limits
    xlim = ax.XLim; % store original xaxis limits
    flag = struct();
    % Transmission values
    % The transmission must be a number between 0 and 1. It is physically
    % impossible for the transmission to be negative or larger than 1 (one).
    
    % BAD FLAG #1
    % If your transmission values are < 0.10(or 10%), the water is too turbid. Disregard these data.
    
    % BAD FLAG #2
    % If transmission shows up as being  larger than 1 (one), then your
    % measurement is most likely taken in very clear water and/or you have a
    % bad zscat measurement obtained with dirty water and/or dirty windows.

    % If your transmission values are > 0.995 (or 99.5%), the water is too clear. Disregard these data.
    
    % QUESTIONABLE FLAG #1
    % Be wary of data collected at transmission values between 0.30 and
    % 0.10 -- generally decreasing data quality as the transmission decreases
    % below 0.30 (30%).
    
    % QUESTIONABLE FLAG #2
    % Be wary of data collected at transmission values between 0.98 and 0.995
    % –low signal-to-noise ratio.
    
    % March 2021 - don't show lower flags
    %flag.values = [0 0.10; 0.995 1.5; 0.10 0.30; 0.98 0.995];
    %flag.name   = {'bad: water too turbid'; 'bad: water too clear'; 'questionable: data quality'; 'questionable: low signal-to-noise'};
    flag.values = [0.995 1.5;  0.98 0.995];
    flag.name   = {'bad: water too clear'; 'questionable: low signal-to-noise'};
   
    for nflag = 1:numel(flag.name)
      % Create temporary variable with flag limits
      tmp = repmat(flag.values(nflag,:),10,1);
      if contains(flag.name{nflag},'bad')
        flag_clr = 'r';%[0 0 0.75];
      else
        flag_clr = 'y';% [0.85 0.85 0.85];
      end
      
      switch which_ax
        case 'y'
          x = linspace(ax.XLim(1),ax.XLim(2),10);
          xfill = [x, x(end:-1:1)]';
          fill(ax,xfill,tmp(:),flag_clr,'FaceAlpha',0.1,'LineStyle','none','DisplayName',flag.name{nflag});
        case 'x'
          y = linspace(ax.YLim(1),ax.YLim(2),10);
          yfill = [y, y(end:-1:1)]';
          fill(ax,tmp(:),yfill,flag_clr,'FaceAlpha',0.1,'LineStyle','none','DisplayName',flag.name{nflag});
      end
    end
    
    % reset axis
    ax.YLim = ylim;
    ax.XLim = xlim;

  end %% FUNCTION PLOT_TRANSMISSION_FLAGS
  close all
end %% MAIN FUNCTION