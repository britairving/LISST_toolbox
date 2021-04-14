function [cfg, data_proc,meta_proc] = LISST_data_qaqc_auto(cfg,data_proc,meta_proc)
%function LISST_data_qaqc_auto runs quality checks on data
%
%  Syntax:
%    [cfg, data_proc] = LISST_data_qaqc_auto(cfg,data_proc)
%
%  Description:
%    Performs automated quality control on processed LISST data. Performs
%    qc tests descripted in manual, as well as basic number could test from
%    in InLineAnaylsis processLISST.m script.
%
%  Refereces:
%    LISST-Deep-Users-Manual-May-2013.pdf "STEP BY STEP PROCEDURE: DATA QUALITY CONTROL"
%    https://seabass.gsfc.nasa.gov/archive/MAINE/boss/EXPORTS/exportsnp/documents/EXPORTS-EXPORTSNP_InLine-LISST-Processing_R1.pdf
%
%  Notes:
%    Incorporate in the future...
%      http://www.sequoiasci.com/article/how-accurate-is-the-concentration-measurement-by-lisst-instruments/
%      https://www.sequoiasci.com/article/the-influence-of-particles-outside-the-size-range-of-the-lisst/
%      https://www.sequoiasci.com/article/interpreting-particle-data/
%    Review and incorporate ambient light effects Andrews et al 2011 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010WR009841
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>

%% 1 | Generate filename for qc'd data
cfg.path.file_qc = strrep(cfg.path.file_proc,'.mat','_qc.mat'); % processed data will be saved here
flag = cfg.qaqc_options.quality_flag; % Pull out flag just for shorthand

%% 2 | Initialize flag array and set all to not evaluated
data_proc.quality_flag        = flag.not_evaluated*ones(size(data_proc.date));
data_proc.quality_flag_test = repmat({''},size(data_proc.date));
% add to metadata
meta_proc.quality_flag.name      = 'data quality flag';
meta_proc.quality_flag_test.name = 'reason for data quality flag, if anything besides not_evalulated or good';

%% 3 | Generate flag metadata 
flags = fieldnames(flag);
meta_proc.quality_flag.flag_meanings = {};
meta_proc.quality_flag.flag_values   = [];
for n = 1:numel(flags)
  meta_proc.quality_flag.flag_meanings = [meta_proc.quality_flag.flag_meanings, flags{n}];
  meta_proc.quality_flag.flag_values   = [meta_proc.quality_flag.flag_values, flag.(flags{n}) ];
end

%% 4 | Loop through automatic QC tests
% QC tests defined in LISST_processing_config.m
% This only works if the tests are set up to be executable with eval and
% fit within the fields of cfg.qaqc_options.quality_flag. 
tests = fieldnames(cfg.qaqc_options.test);
for nqc = 1:numel(tests)
  qctest_name = tests{nqc};
  qctest = cfg.qaqc_options.test.(qctest_name);
  try
    idx_flagged = eval(qctest.method);
  catch
    fprintf(' Could not execute QC test: %s\n',qctest_name)
    fprintf(' **rewrite so it is executable with eval(%s)\n',qctest.method)
  end
  data_proc.quality_flag(idx_flagged)      = flag.(qctest.flag_type);
  data_proc.quality_flag_test(idx_flagged) = {qctest_name};
end

%% 5 | Save qc'd data to .mat file
% if ~ cfg.testing % only save if not in testing mode
fprintf('Saving LISST QC data to %s\n',cfg.path.file_qc)
save(cfg.path.file_qc,'cfg','data_proc','meta_proc');

end %% MAIN FUNCTIOn
