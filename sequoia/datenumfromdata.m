%Function to compute date in a year from columns 39 and 40 from the .dat
%file (read using tt2mat). User must input columns 39 and 40 as an nx2
%matrix, together with the year. Output is the date number in standard
%MATLAB format
%Usage:
%
%   function realdatenum = datenumfromdata(year,data3940);
%
% where:
%   year is the year of the very first measurement
%   data3940 is a nx2 matrix with the first column equal to column 39
%   and the second column equal to column 40 from the .dat file (e.g.
%   data(:,39:40)
%
%OAM, October 30, 2008.

function realdatenum = datenumfromdata(year,data3940)
if nargin~=2
error('There must be 2 input arguments: year and a matrix with 2 columns')
end
if size(data3940,2)~=2;
error('There must be 2 columns in the input matrix data3940')
end

%compute days, hours, minutes, seconds
days = fix(data3940(:,1)/100);
hours = data3940(:,1)-100*days;
minutes = fix(data3940(:,2)/100);
seconds = data3940(:,2)-100*minutes;

years = year*ones(size(data3940,1),1);%creat a vector with the initial year
NewYear=find(diff(days)<0);%find negative differences in daynumber (indicate deployment over new year)
if ~isempty(NewYear);%If we have a new year deployment…
years(NewYear+1:end)=year+1;%…the year after new year is one higher than the inital year (we assume that the deployment doesn’t span more than one new year)!
end

realdatenum = datenum(years,1,days,hours,minutes,seconds);