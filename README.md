# Welcome to the LISST_toolbox!
This toolbox was created to read, process, plot, and write data from profiling LISST (Laser In-Situ Scattering and Transmissometry) instruments manufactured by Sequoia Scientific, Inc.

Currently, the toolbox has only been used to process LISST-Deep Type B and Type C instruments. 

***
# Setup
See Setup wiki page for details on how to organize your data and configuration scripts to utilize this toolbox. 

https://github.com/britairving/LISST_toolbox/wiki/Setup

# Table of Contents
Before running this toolbox, you must organize your data as described in LISST_processing_config.m and configure your project in a script [project]\_config.m, where "project" is the unique name of your project. See the setup wiki page for more details.
## LISST_processing_workflow.m 
wrapper script that calls all subsequent scripts

### [project]\_config.m
"project" is the unique name of your project. For example, LISST_[serial_number]\_[year]\_[experiment]\_[cruiseID] so for the 2018 EXPORTS experiment in the North Pacific, the project name is LISST_sn4025_2018_EXPORTS_SR201812. 
Define metadata used in file header (SeaBASS), instrument specifications, ctd type, and available background files. 

### LISST_processing_config.m 
Builds paths to LISST data and instrument files based on expected folder hierarchy.
Reads/defines instrument information important for data processing and metadata storage.
Aside from path_length, and inversion model used to process data, variables are read from instrument files lisst.ini and InstrumentData.txt using paths built in the beginning part of this script.
Defines data processing model (spherical or random shaped particles), and grid settings.
Defines QC flags and defines automatic QC tests.

### LISST_read_raw_data.m 
Defines structures "meta_raw" and "data_raw".
Converts raw data (*.DAT) files in cfg.path.dir_raw using Sequoia's tt2mat.m script.

### LISST_write_Level1.m
Writes raw data in ascii format.

### LISST_preprocess_data.m
* Defines structures "meta_pre" and "data_pre"
* Calculates date (matlab's datenum format) 
* Calculates temperature and depth based on lisst.ini file contents
* read_ctd_data_by_type.m | Reads in CTD data (subsequent steps much easier if this is upcast and downcast data with true timestamp!)
* LISST_match_ctd_cast.m | Matches LISST profile to CTD cast #
* LISST_correct_time_depth_lag.m | Corrects LISST depth and time based on CTD data
* LISST_identify_downcast.m | Identifies limits of downcast (decending) in LISST profile
* Saves data to cfg.path.file_pre

### LISST_process_data.m
* Defines structures "meta_proc" and "data_proc"
* Limits data to downcasts only
* LISST_select_background_scatterfiles.m | Select which background files to use for processing
* getscat_v2.m | Sequoia's function to process raw data to get scattering and transmission
* invert.p | Sequoia's function that inverts raw data to get uncalibrated volume distribution and midpoint of size bins in microns.
* Calculates other parameters such as volume distribution, particle size distribution, mean particle size, forward scatter coefficient, beam attenuation, etc. 

### LISST_data_qaqc_auto.m
Performs automated quality control on processed LISST data. Performs qc tests described in manual and defined in LISST_processing_config.m, as well as basic number distribution test from in InLineAnaylsis processLISST.m script.

### LISST_grid_data.m
* Defines structures "meta_grid" and "data_grid"
* Grids data (downcasts only) to depth bins defined by cfg.grid_options.bin_depth_m
### LISST_data_qaqc_expert.m
_This is not ready yet but will be done in the future!_

### LISST_write_Level2.m 
* Copies raw files (.dat) and background scatterfiles (.asc) to submit directory
* writes QC'ed data full resolution and gridded data to formatted files based on cfg.write_format
***

# Citation/acknowledgement
If you use this toolbox, please cite as follows:
_Brita Irving, LISST_toolbox, (2021), GitHub repository, https://github.com/britairving/LISST_toolbox/_
***

# References & Resources
* LISST-Deep Users Manual May 2013

* Emmanuel Boss's Ocean Optics package "InLineAnalysis" package (includes flow-through LISST processing from 2018 exports cruise). https://github.com/OceanOptics/InLineAnalysis.git

* Agrawal, Yogesh C., 2005. The optical volume scattering function: Temporal and vertical variability in the water column off the New Jersey coast. Limnology and Oceanography, 6, doi: 10.4319/lo.2005.50.6.1787.

* Andrews, S. W., Nover, D. M., Reuter, J. E., and Schladow, S. G., 2011a. Limitations of laser diffraction for measuring fine particles in oligotrophic systems: Pitfalls and potential solutions. Water Resour. Res., 47, W05523, doi:10.1029/2010WR009837.

* Andrews, S. W., Nover, D. M., Reardon, K. E., Reuter, J. E., and Schladow, S. G., 2011b. The influence of ambient light intensity on in situ laser diffractometers. Water Resour. Res., 47, W06509, doi:10.1029/2010WR009841. 

* Sequoia Scientific, Inc., 2008. Measuring VSF with LISST-100. http://www.sequoiasci.com/article/measuring-absolute-vsf-with-lisst-100/

* Sequoia Scientific, Inc., 2015. Processing LISST-100 and LISST-100X data in MATLAB. https://www.sequoiasci.com/article/processing-lisst-100-and-lisst-100x-data-in-matlab/.
* Jean-Pascal Rueff (2021). Click-n'-drag plot https://www.mathworks.com/matlabcentral/fileexchange/5751-click-n-drag-plot, MATLAB Central File Exchange. Retrieved June 14, 2021.
