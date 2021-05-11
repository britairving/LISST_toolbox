function [cfg, data_grid, meta_grid] = LISST_grid_data(cfg,data_proc,meta_proc)
%FUNCTION LISST_GRID_DATA
%
%  Syntax:
%     [cfg, data_grid] = LISST_grid_data(cfg,data_proc)
%
%  Description:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Set up basic filenames and paths

if cfg.grid_options.ignore_flagged
  cfg.path.file_grid = strrep(cfg.path.file_qc,  '.mat','_gridded.mat'); %fullfile(cfg.path.project,[cfg.project '_processed_qc_gridded.mat']);
else
  cfg.path.file_grid = strrep(cfg.path.file_proc,'.mat','_gridded.mat'); %fullfile(cfg.path.project,[cfg.project '_processed_qc_gridded.mat']);
end

%% 1 | Define number of casts
[unique_casts,u_idx] = unique(data_proc.cast);
if cfg.testing
  num_casts = 6;
else
  num_casts = numel(unique_casts);
end
unique_casts = unique_casts(1:num_casts);

%% 2 | Define which variables to bin
vars_to_bin = ...
  {'date';...
  'lat';...
  'lon';...
  'cast';...
  'station';...
  'botdepth';...
  'depth';...
  'tau';...
  'scat';...
  'cscat';...
  'temp';...
  'VD';...
  'tot_vol_concentration';...
  'mean_size';...
  'PSD';...
  'VSD';...
  'VSF';...
  'PSD_DNSD';...
  'PSD_DVSD';...
  'surf_area';...
  'silt_dens';...
  'silt_vol';...
  'bf';...
  'beam_attenuation';...
  'absorption';...
  'laserPower';...
  'quality_flag'};

% remove variables that don't exist
rm_vars = ~ismember(vars_to_bin,fieldnames(data_proc));
vars_to_bin(rm_vars) = [];
% cast, date, datetime are all gridded differently so remove
% Do not average these per depth bin because will not be monotonic
variables_to_skip = {'datfile' 'zscfile' 'cast' 'lat' 'lon' 'depth' 'station' 'botdepth'};
rm_vars = ismember(vars_to_bin,variables_to_skip);
vars_to_bin(rm_vars) = [];
num_variables = numel(vars_to_bin);

%% 3 | Compute depth intervals
depth_resolution = cfg.grid_options.bin_depth_m;
depth_min   = 0; % round(min(data_proc.depth) / depth_resolution) * depth_resolution;
depth_max   = round(max(data_proc.depth) / depth_resolution) * depth_resolution;
depth_range = depth_min : depth_resolution : depth_max;
num_levels  = numel(depth_range);

%% 5 | Initialize output
data_grid = struct();
data_grid.datfile  = data_proc.datfile(u_idx);
data_grid.zscfile  = data_proc.zscfile(u_idx);
data_grid.cast     = unique_casts(:);
data_grid.lat      = data_proc.lat(u_idx);
data_grid.lon      = data_proc.lon(u_idx);
data_grid.depth    = depth_range;
data_grid.station  = data_proc.station(u_idx);
data_grid.botdepth = data_proc.botdepth(u_idx);

data_grid.bincount = nan(num_casts, num_levels);

% add to metadata
meta_grid.datfile = meta_proc.datfile;
meta_grid.zscfile = meta_proc.zscfile;
meta_grid.cast = meta_proc.cast;
meta_grid.lat  = meta_proc.lat;
meta_grid.lon  = meta_proc.lon;
meta_grid.depth = meta_proc.depth;
meta_grid.depth.grid_resolution = depth_resolution;
meta_grid.depth.grid_min = depth_min;
meta_grid.depth.grid_max = depth_max;
meta_grid.station  = meta_proc.station;
meta_grid.botdepth = meta_proc.botdepth;

meta_grid.bincount.name = 'Number of records averaged into a bin or reported measurement';
meta_grid.bincount.unit = 'none';

for nvar = 1:num_variables
  var = vars_to_bin{nvar};
  size_2 = size(data_proc.(var),2);
  if size_2 == 1
    data_grid.(var) = nan(num_casts, num_levels);
  else
    data_grid.(var) = nan(num_casts, num_levels, size_2);
  end
  % Initialize metadata
  meta_grid.(var) = meta_proc.(var);
  meta_grid.(var).name = ['gridded ' meta_grid.(var).name];
end

%% 6 | Grid data
fprintf('Gridding variables with settings:\n');
fprintf(' depth level min : %d\n', depth_min);
fprintf(' depth level max : %d\n', depth_max);
fprintf(' depth level step: %d\n', depth_resolution);
fprintf(' number of depth levels: %d\n', num_levels);
fprintf(' number of casts       : %d\n', num_casts);
fprintf(' number of variables   : %d\n', num_variables);
rm_casts = [];
fprintf('\n')
for nc = 1:numel(unique_casts)
  ncast = unique_casts(nc);
  fprintf(' Gridding cast %d\n',ncast)
  if cfg.grid_options.ignore_flagged
    cast_select = find(data_proc.cast == ncast & data_proc.quality_flag <= cfg.qaqc_options.quality_flag.not_evaluated);
  else
    cast_select = find(data_proc.cast == ncast);
  end
  if isempty(cast_select) % Speed up when there are no variables
    rm_casts = [rm_casts; ncast];
  else
    cast_depth  = data_proc.depth(cast_select);
    % Loop through variables
    try
      for nvar = 1:num_variables
        var = vars_to_bin{nvar};
        size_2 = size(data_proc.(var),2);
        % get bincount too
        data_grid.bincount(nc,:) =...
          cell2mat(arrayfun(@(d) numel(find(abs(cast_depth-d)<=0.5*depth_resolution)),   depth_range(:), 'UniformOutput', false));
        
        % Nx1 array
        if size_2 == 1
          cast_variables = data_proc.(var)(cast_select,:);
          % also get standard deviation
          if ~strcmp(var,'date')
            data_grid.(var)(nc,:) =...
              cell2mat(arrayfun(@(d) nanmean(cast_variables(abs(cast_depth-d)<=0.5*depth_resolution), 1),   depth_range(:), 'UniformOutput', false));
            data_grid.([var '_sd'])(nc,:,:) =...
              cell2mat(arrayfun(@(d) std(cast_variables(abs(cast_depth-d)<=0.5*depth_resolution, :), 1), depth_range(:), 'UniformOutput', false));
          else
            % smooth top 15 meter so little jogs at top of cast don't cause
            % data_grid time to not monotomically increase
            cast_depth_smo = cast_depth;
            cast_depth_smo(cast_depth_smo <=15) = smooth(cast_depth_smo(cast_depth_smo <=15),15);
            data_grid.(var)(nc,:) =...
              cell2mat(arrayfun(@(d) nanmean(cast_variables(abs(cast_depth_smo-d)<=0.5*depth_resolution), 1),   depth_range(:), 'UniformOutput', false));
          end
        else % % NxM array ( e.g. PSD, VD, VSF, etc )
          cast_variables = data_proc.(var)(cast_select,:,:);
          data_grid.(var)(nc,:,:) =...
            cell2mat(arrayfun(@(d) nanmean(cast_variables(abs(cast_depth-d)<=0.5*depth_resolution, :), 1), depth_range(:), 'UniformOutput', false));
          % also get standard deviation
          for ns = 1:size(cast_variables,2)
            data_grid.([var '_sd'])(nc,:,ns) =...
              cell2mat(arrayfun(@(d) std(cast_variables(abs(cast_depth-d)<=0.5*depth_resolution, ns), 1), depth_range(:), 'UniformOutput', false));
          end % loop through each size bin/ring detector number
        end % if nx1 (e.g. tau) or nxm (e.g. PSD)
      end % loop through variables
    catch
      
    end
  end % if cast has data
  
end % loop through casts

if ~isequal(num_casts,numel(data_grid.cast))
  fprintf('** Expected number of casts / files do not match!\n')
  keyboard
end

%% 8 | Save gridded data to .mat file
fprintf('Saving data_grid data to %s\n',cfg.path.file_grid)
save(cfg.path.file_grid,'cfg','data_grid','meta_grid','data_proc','meta_proc','-v7.3');

end %% MAIN FUNCTION