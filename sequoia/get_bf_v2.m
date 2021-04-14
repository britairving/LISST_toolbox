%  Function to estimate bf directly
%  Created 10/06/03 YCA
%  Mod to include path and feedback resistor type 02/20/04 YCA
%  Edited  Feb 202 Brita Irving to call after getscat already ran
%                  Edited from http://www.sequoiasci.com/article/get_bf-m/
%  Based on paper by YCA in L&O: The optical volume scattering function: temporal and vertical variability inthe water column off the New Jersey Coast. Limnology and Oceanography 50: 1787–1794.
%  http://www.sequoiasci.com../graphics/upload/pdf/library/technicalpapers/Agrawal_2005_LO.pdf
%
%  Usage:
%    [bf,c,a,laserPowerEnteringWater]=get_bf_v2(tau,data,cscat,ringareafile,HK3Scale,HK4Scale,HK4Off,path,X);
%
%   where
%   datafile:         Your LISST-100/-100X .DAT datafile
%   zscfile:          Your LISST-100/-100X zscat file, usually with .ASC extension (if obtained with the LISST-SOP software)
%   dcal_file:        Your LISST-100X RingArea_xxxx.ASC file, where xxxx is serial number
%   HK3Scale:         mW laser power entering water per digital count in factory zscat variable 36; Look it up in LISST.INI for your serial #.
%   HK4Scale,HK4Off:  pressure scales and offsets for calibration of pressure sensor. Look it up in LISST.INI for your serial #.
%   path_cm:          path length in cm, usually 5 for a LISST-100/-100X
%   X:                set X = 1 for a LISST-100X, any other value for a LISST-100
%
%   OAM edits 12/09/08 to correct error for dividing by zscat laser reference when computing P0.
%   OAM edits 4/4/11 to minor unnecessary pre-allocations etc.
%

function [bf,c,a,laserPowerEnteringWater]=get_bf_v2(tau,data,cscat,ringareafile,HK3Scale,HK4Scale,HK4Off,path,X)
[r,sc]=size(cscat);
dcal=load(ringareafile);%load ringarea file
if(r==40 && sc==1)%if there’s only one measurement
  dcal=dcal';
end

laserPowerEnteringWater = data(:,36).*HK3Scale;%this is now laser power entering water, in units of mW

press = data(:,37)*HK4Scale+HK4Off; %compute the depth in calibrated units of m using factory supplied constants from LISST.INI file
mins = fix(data(:,40)/100); secs=data(:,40)-100*mins; days=fix(data(:,39)/100); hours=data(:,39)-100*days; %compute time
time = days+hours/24+mins/(24*60)+secs/(24*60*60); %create decimal date/time

%pre-allocate
bf=ones(r,1);

% Compute forward scattering coefficient up to max range covered by LISST.
for i=1:r
  bf(i)=3.2e-4*sum(cscat(i,:)); % division by tau already done in computing scat; at this stage b has units of mW/m. Constant of 3.2e-4 is from combining detector gains with resistivity and silicon sensitivity
  bf(i)=(bf(i)/laserPowerEnteringWater(i))*5/path;% units are 1/m after dividing by P0. It is assumed that the path length is 5 cm, so adjust wiht proper ratio if this is not the case.
end
%compute beam attenuation, again adjusting for cases where the path is different from 5 cm.
% 1/0.05m -->20m-1
% SeaBASS c [1/m]  = Beam attenuation coefficient 
% SeaBASS cp [1/m] = Beam attenuation coefficient of particles (ap + bp)
% <https://seabass.gsfc.nasa.gov/wiki/stdfields#Table%20of%20Field%20Names%20and%20Units>
c=-20*log(tau).*5/path;
%3.2e-4 assumes 5V/4095 counts – for LISST-100X it is different so correct
%(if necessary) as follows
if X==1 %In case of LISST-100X
  volts = 2.5;
else %in case of LISST-100
  volts = 5;
end
bf=bf*volts/5;% Correct if necessary
a=c-bf;%In principle this is now absorption