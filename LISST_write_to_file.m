function LISST_write_to_file(cfg, data_proc, meta_proc, data_grid, meta_grid)
%FUNCTION LISST_write_to_file
%
%  Syntax:
%    [cfg, dat] = LISST_write_to_file(cfg,dat)
%  
%  Description:
%
%  Refereces:
%     https://seabass.gsfc.nasa.gov/wiki/Data_Submission
%     https://seabass.gsfc.nasa.gov/wiki/stdfields#Table%20of%20Field%20Names%20and%20Units
%
%  Notes:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
% datfile = cfg.datfile{1};
% if ismac && contains(datfile,'\')
%   rm_str = strfind(datfile,'\');
%   winpath = datfile(1:rm_str(end));
%   cfg.datfile = erase(cfg.datfile,winpath);
% end
% %% 1 | Create submit folder with specifics on inversion model and background files
% if ~isfield(cfg.path,'submit_dir')
%   if cfg.inst.rand
%     cfg.path.submit_dir = fullfile(cfg.path.project,['submit_nonspherical_' cfg.zscat_choice]);
%   else
%     cfg.path.submit_dir = fullfile(cfg.path.project,['submit_spherical_' cfg.zscat_choice]);
%   end
%   cfg.path.submit_doc = fullfile(cfg.path.submit_dir,'documents'); 
%   mkdir(cfg.path.submit_dir) % store all formated data files
%   mkdir(cfg.path.submit_doc) % store associated documents/zscat files
% end

% %% 2 | move matlab datafile, zscat files, and raw data files into submit/documents folder
% fprintf('Copying raw data files and background files to %s\n',cfg.path.submit_doc)
% copyfile(cfg.path.fname_binned,cfg.path.submit_doc)
% for ndatf = 1:numel(cfg.datfile)
%   if ~exist(fullfile(cfg.path.submit_doc,cfg.datfile{ndatf}),'file')
%     copyfile(fullfile(cfg.path.raw,cfg.datfile{ndatf}),cfg.path.submit_doc);
%   end
% end
% 
% unique_zscfiles = unique(cfg.zscfile);
% for nzscf = 1:numel(unique_zscfiles)
%   if ~exist(fullfile(cfg.path.submit_doc,unique_zscfiles{nzscf}),'file')
%     copyfile(fullfile(cfg.path.zscdir,unique_zscfiles{nzscf}),cfg.path.submit_doc);
%   end
% end


switch cfg.write_format
  %% I  | SeaBASS format
  case 'seabass' % https://seabass.gsfc.nasa.gov/wiki/Data_Submission
     LISST_write_SeaBASS(cfg,data_proc,data_grid,meta_proc)
  %% II | NGA LTER Format
  case 'nga_lter' % AXIOM Research Workspace
    LISST_write_NGALTER(cfg,dat,dat_bin,meta)
  case 'netcdf'
    fprintf('This write format is set up yet\n')
    keyboard
  otherwise
    fprintf('This write format is set up yet\n')
    keyboard
end

end %% MAIN FUNCTION 