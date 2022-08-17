function zscat = LISST_select_insitu_background(cfg,downcasts,zscat)
%function LISST_select_background_scatterfiles
%
%  Syntax:
%    cfg = LISST_select_background_scatterfiles(cfg)
%
%  Description:
%
%  Refereces:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% TO DO
fprintf('INCORPORATE NEW INFO FROM SEQUOIA!!\n')
%https://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/
% III.2  Using a local data record as a background file
% 
% This is often useful in particularly sensitive environments. For example,
% if you are doing a profile in the water column and suspect that the
% instrument is going out of alignment due to pressure (we admit, sometimes
% this may happen!), then find a local minimum again using the ginput
% command, or by plotting some other variable. Example methods may be: find
% the point where t is highest. This can be done by:

%     [a,b] = max(tau);
%     zsc = data(b,:); 
% 
% One may then save this estimate of the background and process the data as
% already described. To save, type
% 
%     save zsc.asc zsc ?ascii 
% 
% Another example is encountered in looking at deep ocean profiles. Here,
% instead of the maximum transmission, one may find a local point that is
% suitable as a background. For example, plot the ratio
% data(:,33)/data(:,36) against depth [data(:,37)]. This is something like
% optical transmission. The data may involve profiles across a hydrothermal
% plume. Find a point just above the plume and call it the local
% background. You may use the ginput command to identify this local
% background. Again, you may save this as a background, or use it in Matlab
% in the detailed procedure described in II.2.
%%
use_method = 2;  % 
z_window   = 10; % m
%% Try using an “in-situ background”
% The fact that the in-situ transmission goes slightly above
% 1 is due to normal small imperfections in the measurement, and indicates
% that the in-situ water is virtually pure, at least as far as the LISST
% can tell. Unfortunately, that also means there were simply too few
% particles for the LISST to make useful size measurements. I would not
% trust particle size distributions derived from this dataset, at least not
% with standard processing. 
% But you could try using an ?in-situ background?, that is, find a section
% of data that has consistently higher transmission than the rest, average
% it, and then use that as your background. Sometimes that will yield
% useful PSDs for the segments of data that have lower transmission. But if
% you do this, you will need to look skeptically at the results and
% consider whether they are physically reasonable.

%% Plot transmission for each profile
% transmissions calculated with factory background
ax = makefig_subplots(2,1);
plot(ax(1),zscat.factory.tau.time_m,zscat.factory.tau.tau_m,'.');
hold(ax(1),'on'); grid(ax(1),'on');
plot(ax(1),zscat.factory.tau.time_m,nanmean(zscat.factory.tau.tau_m,2),'k-','LineWidth',4);
datetick(ax(1),'x','keeplimits','keepticks')

plot(ax(2),zscat.factory.tau.tau_m,zscat.factory.tau.depth_m);
ax(2).YDir = 'rev'; grid(ax(2),'on'); hold(ax(2),'on');

%ax(1).YLim = [0.85 0.95];
%ax(2).XLim = ax(1).YLim;
ax(2).YAxisLocation = 'right';
ax(1).XTickLabelRotation = 45;

ax(2).XLabel.String = 'Transmission';
ax(1).XLabel.String = 'Date';
ax(1).YLabel.String = 'Transmission';
ax(1).Title.String  = strrep(cfg.project,'_','\_');
ax(2).Title.String  = 'In-situ background cast selection';

%% Find the cast with highest transmission values 
% There are multiple ways to do this... adding another subjective component
% to the data processing!

% For now, just select the profile with the highest average transmission
% between a depth range, then use the highest transmission point to
% average and use as the background
if strcmp(cfg.project,'LISST_sn4041_2021_EXPORTS_JC214') || strcmp(cfg.project,'LISST_sn4025_2021_EXPORTS_DY131')
  depth_range = [1500 2000];
else
  depth_range = [100 500];
end

zrng = find(zscat.factory.tau.depth_m > depth_range(1) & zscat.factory.tau.depth_m < depth_range(2));
switch use_method
  case 1
    % Method 1: find cast with highest average transmission between depth ranges
    tau_mean = nanmean(zscat.factory.tau.tau_m(:,zrng),2);
    [~,imax_mean] = max(tau_mean,[],1);
    plot(ax(2),zscat.factory.tau.tau_m(imax_mean,:),zscat.factory.tau.depth_m,'g-','linewidth',4);
    castnum = zscat.factory.tau.cast_m(imax_mean);
    [~,datname,~] = fileparts(cfg.proc_options.datfile{imax_mean});
    % find depth of highest transmission 
    [~,idx_z_maxtau] = max(zscat.factory.tau.tau_m(imax_mean,:));
    depth_at_tau_max = zscat.factory.tau.depth_m(idx_z_maxtau);
    fprintf('Maximum average tau values in %s\n',datname)
    fprintf('Maximum tau value at %d m\n',depth_at_tau_max)
  case 2  
    % find single highest transmission point
    [max1,~] = max(zscat.factory.tau.tau_m(:,zscat.factory.tau.depth_m > depth_range(1)),[],2);
    [~,imax2] = max(max1);
    plot(ax(2),zscat.factory.tau.tau_m(imax2,:),zscat.factory.tau.depth_m,'k-','linewidth',5);
    castnum = zscat.factory.tau.cast_m(imax2);
    % Method 2: find cast with highest transmission above maximum depth
    [~,izmax] = max(zscat.factory.tau.tau_m(imax2,zrng));
    idx_z_maxtau = zrng(izmax);
    depth_at_tau_max = zscat.factory.tau.depth_m(idx_z_maxtau);
    plot(ax(2),zscat.factory.tau.tau_m(imax2,idx_z_maxtau),zscat.factory.tau.depth_m(idx_z_maxtau),'ks','markerfacecolor','g','markersize',10);
    
    [~,datname,~] = fileparts(cfg.proc_options.datfile{imax2});
    fprintf('Maximum average tau values in %s\n',datname)
    fprintf('Maximum tau value at %d m\n',depth_at_tau_max)
  otherwise
    fprintf('Not set up yet...\n')
    keyboard
end

if cfg.savefig
  figname = fullfile(cfg.path.dir_figs,[cfg.project '_background_insitu_cast_selection']);
  standard_printfig_lowrespng(figname)
end
%% Update depth_range to be within +/- 10m around highest transmission
depth_range = [depth_at_tau_max-z_window depth_at_tau_max+z_window];

%% smooth data to get background around depth_range
% Fix issue with paths 
if ismac && contains(datname,'\')
  rm_str = strfind(datname,'\');
  datname = datname(rm_str(end)+1:end);
end

% Limit to downcast and between depths
pidx = find(downcasts.cast == castnum);
depth = downcasts.depth(pidx);
pidx2 = find(depth > depth_range(1) & depth < depth_range(2));
% Reindex to original resolution
pidx = pidx(pidx2);

% Pull out unique depth values so can interpolate
[dep,iux] = unique(downcasts.depth(pidx));
zbin = floor(dep(1)):5:ceil(dep(end));
data = nan(40,numel(zbin));
try
  for ii = 1:40
    tmp = downcasts.data(pidx(iux),ii);
    tmp = smoothdata(tmp,'movmean',17);
    data(ii,:) = interp1(dep,tmp,zbin,'linear','extrap');
    %data(ii,:) = smoothdata(data(ii,:),'movmean',5);
  end
catch
  keyboard
end

%% Add to zscat structure
zscat.insitu.file = [datname '_zsc_average' num2str(depth_range(1)) '_to_' num2str(depth_range(2)) 'm.asc'];
zscat.insitu.date = nanmean(downcasts.date(pidx(iux)));
zscat.insitu.data = nanmean(data,2);
% Calculate r = the laser power/laser reference ratio (to adjust for drift in laser output power over time)
r = zscat.insitu.data(33)/zscat.insitu.data(36);
zscat.insitu.tau = LISST_calculate_transmission(cfg,downcasts,r);

insitu_zscat_file = fullfile(cfg.path.dir_zsc,zscat.insitu.file);
fprintf('Saving insitu background to %s\n', insitu_zscat_file)
fileID = fopen(insitu_zscat_file,'w');
fprintf(fileID,'%0.4f\n',zscat.insitu.data);
fclose(fileID);


end %% MAIN FUNCTION