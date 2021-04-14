function ax = makefig_subplots(subplotsx,subplotsy)
% make standard figure
% taken from http://p-martineau.com/perfect-subplot-in-matlab/

%parameters for figure and panel size
% make new figure
try 
  uname = char(java.lang.System.getProperty('user.name'));
catch
  
end
if exist('uname','var')
  switch uname
    case 'Rachel'
      plotheight = 22;
      plotwidth  = 17;
      pos_x      = 0.80;
      pos_y      = 9.5;
    case 'bkirving'
      plotheight = 22;
      plotwidth  = 17;
      pos_x      = 0.80;
      pos_y      = 9.5;
     case 'Stephanie'
      plotheight = 22;
      plotwidth  = 17;
      pos_x      = 0.80;
      pos_y      = 9.5; 
  end
else
  plotheight = 22;
  plotwidth  = 17;
  pos_x      = 0.80;
  pos_y      = 9.5;
end
if nargin == 0
subplotsx = 3;
subplotsy = 5;
end
leftedge=1.2;
rightedge=0.6;   
topedge=1.5;
bottomedge=1.5;
spacex=0.4;
spacey=0.3;
fontsize=14;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
%setting the Matlab figure
% make new figure
n = 1;

while ishandle(n)
  n = n + 1;
end
fig = figure(n);
fig.Units             = 'inches';
fig.Position          = [plotheight pos_x plotwidth pos_y];
fig.PaperPositionMode = 'auto';
fig.PaperOrientation  = 'portrait';
fig.PaperSize         =  [plotwidth plotheight];
fig.PaperPosition     = [0 0 plotwidth plotheight];
 
%loop to create axes
for i=1:subplotsx
  for ii=1:subplotsy
    ax(i,ii)=axes('position',sub_pos{i,ii},'Box','on','Layer','top');%,'XGrid','off','XMinorGrid','off','FontSize',fontsize,);
  end
end

set(fig, 'color', 'w')
set(fig,'DefaultAxesLineStyleOrder','-|--|-.|:');
tomLnCmap = [ 0 0 1; 1 0 0; 1 0.6 0; 0.85 0 0.95; 0 .95 0.95; 0.8 0.8 0; 0.8 0.5 0; 0.5 0 0; 0 0 0.5; 0 0.65 0.65];
%               b;     r;   orange;      ~m;        ~bc;          ~y;    dorange;      dr;     db;       dc;
% no greens as Tom is R-G color blind.
set(fig,'DefaultAxesColorOrder',tomLnCmap);
% axes defaults
set(fig,'DefaultTextFontname',    'Times');
set(fig,'defaultTextFontSize',    16);
set(fig,'DefaultTextFontName',    'Times');
set(fig,'DefaultAxesFontName',    'Times');
set(fig,'DefaultAxesFontSize',    16);
set(fig,'DefaultAxesFontWeight',  'Bold');
end