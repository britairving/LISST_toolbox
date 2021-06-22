function EXPORTS_Epoch_comparison_LISST_sb

SR_sb = '/Volumes/toshiba1/LISST_Data_Structured/LISST_sn4025_2018_EXPORTS_SR201812/submit/L2/EXPORTS-EXPORTSNP_LISST-Deep_survey_20180814-20180909_R0.sb';
[SR, sbHeader, headerArray] = readsb(SR_sb, 'MakeStructure',true);

RR_sb = '/Volumes/toshiba1/LISST_Data_Structured/LISST_sn4041_2018_EXPORTS_RR201813/submit/L2/EXPORTS-EXPORTSNP_LISST-Deep_process_20180814-20180908_R0.sb';
[RR, sbHeader, headerArray] = readsb(RR_sb, 'MakeStructure',true);



end