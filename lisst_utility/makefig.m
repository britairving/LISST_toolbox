function fig = makefig
%% function fig = makefig
% make uvp standard figure 
% author: Brita Irving <bkirving@alaska.edu>
%%
% uses and resets current figure
% fig = figure(1);
% clf('reset');
% make new figure
try 
  uname = char(java.lang.System.getProperty('user.name'));
catch
end
axiscolor = [0.9 0.9 0.9];
n = 1;
while ishandle(n)
  n = n + 1;
end
fig = figure(n);
fig.Units = 'inches';
if exist('uname','var') && strcmp(uname,'bkirving')
  fig.Position = [21 0.5104 17 9.4479];
elseif exist('uname','var') && strcmp(uname,'Brita Irving')
  fig.Position = [11.0833    0.4271    8.9271   10.5000];
else
  fig.Position = [1.3021    0.5104   14.0312    7.2479];
end
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'portrait';
% axes defaults
set(fig,'DefaultTextFontname',    'Times');
set(fig,'defaultTextFontSize',    16);
set(fig,'DefaultTextFontName',    'Times');
set(fig,'DefaultAxesFontName',    'Times');
set(fig,'DefaultAxesFontSize',    20);
set(fig,'DefaultAxesFontWeight',  'Bold');
if ismac
  set(fig,'DefaultFigureRenderer', 'opengl');
end
set(fig, 'color', 'w');
set(fig,'DefaultAxesLineStyleOrder','-|--|-.|:');

tomLnCmap = [ 0 0 1; 1 0 0; 1 0.6 0; 0.85 0 0.95; 0 .95 0.95; 0.8 0.8 0; 0.8 0.5 0; 0.5 0 0; 0 0 0.5; 0 0.65 0.65];
%               b;     r;   orange;      ~m;        ~bc;          ~y;    dorange;      dr;     db;       dc;
% no greens as Tom is R-G color blind.
set(fig,'DefaultAxesColorOrder',tomLnCmap);
% set axis color
% set(gca,'color',[0.9 0.9 0.9])
set(gca,'YDir','reverse','XAxisLocation','bottom','Layer','top','Color',axiscolor);
end