function EXPORTSNA_final_LISST_tidy
% Description: Load processed/gridded file, update any pertinent
% information to cfg so SeaBASS header is correct, then write SeaBASS
% files, and write to csv for easy access via Google Drive
%
% Author: Brita Irving <bkirving@alaska.edu>
%% DY131
DY131 = 'F:\LISST_Data_Structured\LISST_sn4025_2021_EXPORTS_DY131\proc\LISST_sn4025_2021_EXPORTS_DY131_proc_spherical_gridded.mat';
load(DY131);
% update cruise
cfg.header.cruise = 'EXPORTSNA';
% update operator 
cfg.header.comments = strrep(cfg.header.comments,'LISST operator: Jessica Pretty','LISST operator: Rachel Lekanoff');
% update "calibration file" i.e. background file
cfg.header.calfiles = unique(data_grid.zscfile);

% LISST_write_Level2(cfg, data_proc, meta_proc, data_grid, meta_grid);
clear
%% JC214
JC214 = 'F:\LISST_Data_Structured\LISST_sn4041_2021_EXPORTS_JC214\proc\LISST_sn4041_2021_EXPORTS_JC214_proc_spherical_gridded.mat';
load(JC214);
% update cruise
cfg.header.cruise = 'EXPORTSNA';
% update operator 
cfg.header.comments = strrep(cfg.header.comments,'LISST operator: Jessica Pretty, University of Alaska Fairbanks','LISST operator: Jordan Snyder');
% update "calibration file" i.e. background file
cfg.header.calfiles = unique(data_grid.zscfile);
% LISST_write_Level2(cfg, data_proc, meta_proc, data_grid, meta_grid);
clear
