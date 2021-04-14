function [profile_index, profile_direction] = bi_find_downcast(depth_var, time, profile_limit, plotit)
%% ------------------------------------------------------------------------
% function [profile_index profile_direction] = ...
%                              bi_find_downcast(depth_var,time, profile_limit, plotit)
% Description:
%   Automatically detect where full descending profile begins and ends.
%   
%   profile_direction = -1 ; % descending
%   profile_direciotn = +1 ; % ascending
% 
% BKI 02.2017, 09.18.2017
%       profile_limit: minimum length of a profile, profiles are discarded
%         if shorter than profile_limit [dbar].
%       profile_soak_depth
%% ------------------------------------------------------------------------
dbstop if error
if nargin < 4
    plotit = 0;
end
% set the size of time bin to find the turning points
dt = 5./24./60./60; % 5-second bin size
p = fillInvalidValues(depth_var,'linear');
time = time;

profile_index     = nan(size(time));
profile_direction = nan(size(time));
% interpolate over any nans at the ends of the time vector if they exist
time = interp1( find( ~isnan(time)), time(~isnan(time)), 1:length(time), 'linear', 'extrap')';
% do the same thing with pressure
p = interp1( find( ~isnan(p)), p(~isnan(p)), 1:length(p), 'linear', 'extrap')';
% now low pass filter pressure to try and avoid erroneous data
lowpass_span = 11;

p = smooth(p,lowpass_span);

% bin the pressure
tbin = [nanmin(time):dt:nanmax(time)]'; 
pbin = binaverage( time, p, tbin ); 
% interpolate nans
pbin = fillInvalidValues(pbin,'linear');

% use the binned data to calculate first and second derivatives of pressure
dpdt = gradient( pbin, tbin );   % vertical velocity
d2pdt2 = gradient( dpdt, tbin ); % double derivative of pressure

% make an up-down vector from the binned pressure data
% (+1 rising, -1 sinking, 0 constant depth)
udvec = nan*pbin;
udvec( sign( dpdt ) <0 ) = +1; % find all the ascending data points
udvec( sign( dpdt ) >0 ) = -1; % find all the descending data points
udvec(isnan(udvec)) = 0;

% find all the turning points in the up-down vector
tpb = find( diff( udvec) ~=0 ); 

% catch if only one turning point
if numel(tpb) == 1
  tpb = [ tpb; numel(udvec)];
end
% catch if first profile identified as ascending
if tpb > 20
  if udvec(tpb(2)) == +1 && nanmedian(udvec(1:tpb(1))) == -1
    if find(pbin > 0.1,1) < tpb(1)
      tpb = [find(pbin > 0.1,1); tpb];
    end
  end
end

% shift one indice to the right
tpb = tpb+1; 
% now to only keep the turning points that define the start of a new
% profile. Empirically decided that profiles with less than 5 data points
% are small 'jogs' undergone by the glider rather than a total profile.
faketurns = (diff( tpb) < 11); % finplotd fake turns dur to jogs
% marks profile end as just before upcast instead of end of downcast
% to correct for this, shift faketurns by 1
faketurns = circshift(faketurns,1);
tpb(faketurns) = nan;  % throw out fake turns.
tpb = tpb(~isnan( tpb )); 

% Define whether each profile is up or down
updowns = udvec( tpb ); 
% flag bad indices where profiles are same direction 
badinds = find( diff(updowns) ==0 ); 
badinds = badinds+1;
% Make sure not removing turning point before true profile begins.. this
% could happen if long soak or something
[~,imax] = max(pbin(tpb));
before_descent = imax-1;
if ismember(before_descent,badinds)
  badinds(badinds == before_descent) = [];
end
% now remove the bad indices
updowns(badinds) = nan; 
updowns = updowns( ~isnan( updowns )); 
tpb(badinds) = nan; 
tpb = tpb( ~isnan( tpb )); 
% keyboard
%  now use find the turning points in the original CTD data vector

tpCTD = interp1( time, 1:length( time ), tbin(tpb) ); 
tpCTD = floor( tpCTD ); 
tpCTD = [tpCTD; length( time )]; 

% make sure turning points are on finite pressure values
finite_indicies = find(isfinite(depth_var));
if any(isnan(depth_var(tpCTD)))
  fprintf('test this\n')
  keyboard
  for ntp = 1:numel(tpCTD)
    idx_near = nearest(finite_indicies,tpCTD(ntp));
    tpCTD(ntp) = finite_indicies(idx_near);
  end
end

%% Make sure profile_limit is appropriate for cast
if max(p(tpCTD)) < profile_limit
  % find soak depth
  turning_points = p(tpCTD);
  [~,imax] = max(turning_points);
  max_turnpoints = turning_points(imax);
  turning_points(imax) = [];
  soak_depth = max(turning_points);
  % If the longest profile is at least half of the user defined profile
  % limit, then reset the profile limit to catch unusually shallow profiles
  if max_turnpoints - soak_depth > profile_limit/2
    profile_limit = (max_turnpoints - soak_depth) - 1;
  end
end

%% discard profiles that are less than profile_limit [dbar] 
profiles     = length( tpCTD )-1;
badprofiles = [];
for cc = 1:profiles
  % get indicies from turning point to turning point
  rng = tpCTD(cc):tpCTD(cc+1)-1;
  % make sure not just single point
  if length(rng) <= 1
    badprofiles = [badprofiles cc];
    continue
  end
  % make sure profile is greater than profile_limit [dbar]
  
  dbar = abs(p(rng(1)) - p(rng(end)));
  if dbar <= profile_limit
    badprofiles = [badprofiles cc];
    continue
  end
end
profiles = profiles - length(badprofiles);
tpCTD(badprofiles)   = [];
updowns(badprofiles) = [];

%% Fill profile_index and profile_direction 
for np = 1:profiles
  % get indicies from turning point to turning point
  rng = tpCTD(np):tpCTD(np+1)-1;
  profile_index(rng) = np;              % profile identifier
  profile_direction(rng) = updowns(np); % profile direction 
end

%% Want to keep the same format as findProfiles_socib.m
% Therefore, mark all inflection points and surface points as 0
% profile_direction and profile_index with +0.5 from last profile_index
profile_direction(tpCTD) = NaN;
profile_index(tpCTD)     = NaN;

% Now fill profile_direction(isnan) = 0;
profile_direction(isnan(profile_direction)) = 0;
% loop through all indices and fill profile_index( profile_direction == 0)
% with last profile_index + 0.5

ifin = find(isfinite(profile_index));
if isempty(ifin)
  fprintf('No profile identified\n')
  return
end

profile_index(1:ifin(1)-1) = 0.5;
ilast = numel(profile_index);
for i = 1:ilast
  if isnan(profile_index(i))
    bkwdrng = 1:i;
    frwdrng = i:ilast;
    lastidx = find(isfinite(profile_index(bkwdrng)),1,'last');
    lastidx = bkwdrng(lastidx);
    nextidx = find(isfinite(profile_index(frwdrng)),1,'first');
    nextidx = frwdrng(nextidx);
    if isempty(nextidx) % make sure last indices are filled as well
      profile_index( lastidx+1 : ilast ) = profile_index(lastidx) + 0.5;
      continue
    else % fill points between profiles with last profile_index + 0.5
      profile_index( lastidx+1 : nextidx-1 ) = profile_index(lastidx) + 0.5;
    end
  end
end

%% Since using this script to find a single profile in LISST data.. simple catch here in case bump mid cast
% p1 = find(profile_index == 1);
% if max(depth_var(p1)) < max(depth_var)
%   [~,idx_max] = max(depth_var);
%   profile_index(p1(1):idx_max) = 1;
%   profile_direction(p1(1):idx_max) = profile_direction(p1(1));
% end

% by default don't plot the data
if plotit
    % plot it up
    figure1=figure('Position', [1950, 50, 1800, 900]); clf
    plot( tbin, pbin, 'k-'); hold on
    axis ij; datetick
    % plot( tbin(tpb), pbin(tpb), 'kp', 'markerfacecolor', 'r', 'markersize', 10)
    for i = 1:profiles
      idx_prof = find(profile_index == i);
      if isempty(idx_prof)
        continue
      end
      if profile_direction(idx_prof(1)) == -1 % downcast
        plot( time(idx_prof), p(idx_prof), 'b.', 'markersize', 7)
      elseif profile_direction(idx_prof(1)) == +1  % upcast
        plot( time(idx_prof), p(idx_prof), 'r.', 'markersize', 7)
      end
      text(nanmean(time(idx_prof)),nanmean(p(idx_prof)),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','center','BackgroundColor','w','margin',1,'clipping','on')
    end
    titleout( 'final profiles')
    plot( time(tpCTD), p(tpCTD), 'kp', 'markerfacecolor', 'c', 'markersize', 10)
    datetick('x','mm/dd HH:MM','keeplimits')
    fprintf('paused here... press any key to continue\n')
    pause
    close
end
%% make sure first profile is decending
pdown = find(profile_index == 1);
if profile_direction(pdown(1)) == 1
  fprintf('First profile is ascending...\n')
  profile_index(:) = NaN;
  profile_direction(:) = NaN;
  % select manually
  figure('Position', [1950, 50, 1800, 900]); clf
  axes(gca);
  done = 0;
  while ~done
    plot( tbin, pbin, 'k-'); hold on; axis ij; datetick
    title('Manually select profile limits')
    fprintf('Select profile starting place on plot\n')
    plot_select = ginput(1);
    startind = find(time >= plot_select(1),1);
    plot( time(startind), p(startind), 'kp', 'markerfacecolor', 'g', 'markersize', 15)
    fprintf('Select profile ending place on plot\n')
    plot_select = ginput(1);
    endind = find(time <= plot_select(1),1,'last');
    plot( time(endind), p(endind), 'kp', 'markerfacecolor', 'r', 'markersize', 15)
    profile_index(startind:endind) = 1;
    profile_direction(startind:endind) = -1;
    plot( time(startind:endind), p(startind:endind), 'b-', 'LineWidth', 2)
    chc = input('  Is this correct? <1/0> ');
    if chc == 0
      done = 0;
      clf
    else
      close(gcf);
      done = 1;
    end
  end
end
