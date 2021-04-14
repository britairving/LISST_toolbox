function tau = LISST_calculate_transmission(cfg,downcasts,r)
%% FUNCTION LISST_CALCULATE_TRANSMISSION
%
%  Syntax:
%
%  Description:

%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
% initialize tau
tau = struct();
tau.tau   = [];
tau.time  = [];
tau.depth = []; %
tau.cast  = []; % sequential number

%% Pull out unique casts
casts = unique(downcasts.cast,'stable');

%% Calculate transmission
% compute optical transmission, taking the eventual drift in laser power into account
tau.tau   = downcasts.data(:,33) ./r ./ downcasts.data(:,36);
% Pull out time (matlab datenum format)
tau.time  = downcasts.date;
% Pull out depth
tau.depth = downcasts.depth;
% Store the sequential number, not necessarily the actual cast number
tau.cast  = downcasts.cast;


%% Grid and smooth each tau profile
% Create tau matrix of [profiles x depths]
zmax = fix(max(tau.depth));
if zmax > 1000
  zbin = 1:25:zmax;
elseif zmax > 200
  zbin = 1:10:zmax;
else % shallow shelf
  zbin = 1:5:zmax;
end

tau_m   = nan(numel(cfg.proc_options.datfile),numel(zbin));
times   = nan(numel(cfg.proc_options.datfile),1);
castm   = nan(numel(cfg.proc_options.datfile),1);

for nf = 1:numel(casts)
  cast_idx = find(tau.cast == casts(nf));
  times(nf) = tau.time(cast_idx(1));
  castm(nf) = tau.cast(cast_idx(1));
  tmp = smoothdata(tau.tau(cast_idx),'movmean',17);
  % Pull out unique depths, otherwise interp1 will not work
  [dep,iux] = unique(tau.depth(cast_idx)); tmp = tmp(iux);
  % interpolate profile tau to binned depth
  tau_m(nf,:) = interp1(dep,tmp,zbin,'linear');
  % Smooth more
  tau_m(nf,:) = smoothdata(tau_m(nf,:),'movmean',5);
  % get rid of erroneous data
  rm_depths = find(zbin < dep(1) | zbin > dep(end));
  tau_m(nf,rm_depths) = nan;
  % delete last few entries in case false spike
  ifin = find(isfinite(tau_m(nf,:)));
  if numel(ifin) > 10
    tau_m(nf,ifin(end-5:end)) = nan;
  end
end

tau.time_m  = times;
tau.cast_m  = castm;
tau.depth_m = zbin;
tau.tau_m   = tau_m;
% Calculate the average tau profile with depth
tau.mean = nanmean(tau_m,1);
% Calculate the average tau over the top 200m
idx_surf = zbin >= 0 & zbin <= 200;
tau.surf = nanmean(tau_m(:,idx_surf),2);

% Remove the ungridded data, not really necessary at this point
tau = rmfield(tau,'cast');
tau = rmfield(tau,'tau');
tau = rmfield(tau,'time');
tau = rmfield(tau,'depth');

end %% FUNCTION CALCULATE_TRANSMISSION
