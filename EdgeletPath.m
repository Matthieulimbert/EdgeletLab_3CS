% EdgeletPath:

restoredefaultpath;

global EDGELETLABPATH
global PATHNAMESEPARATOR
global PREFERIMAGEGRAPHICS
global MATLABPATHSEPARATOR
		
PREFERIMAGEGRAPHICS = 1;
Friend = computer;

if strcmp(Friend,'MAC2'),
  PATHNAMESEPARATOR = ':';
  EDGELETLABPATH = [pwd, PATHNAMESEPARATOR];
  MATLABPATHSEPARATOR = ';';
elseif isunix,
  % Mac OS X returns isunix=1
  PATHNAMESEPARATOR = '/';
  EDGELETLABPATH = [pwd, PATHNAMESEPARATOR];
  MATLABPATHSEPARATOR = ':';
elseif strcmp(Friend(1:2),'PC');
  PATHNAMESEPARATOR = '\';	  
  EDGELETLABPATH = [pwd, PATHNAMESEPARATOR];  
  MATLABPATHSEPARATOR = ';';
end

disp('-- STARTING EDGELETLAB --')
disp('Adding EdgeletLab to MATLAB path...')

post = PATHNAMESEPARATOR;
p = path;
pref = [MATLABPATHSEPARATOR EDGELETLABPATH];
p = [p pref];

p = [p pref 'EdgeletTrans' post];

p = [p pref 'Networks' post ];

p = [p pref 'Experiences' post];

p = [p pref 'Utilities' post];

p = [p pref 'Images' post];

% mex code
p = [p pref 'mex' post 'src' post 'Networks' post];

path(p);


% Check if MEX files have been compiled and compile them if needed
CompileMex;

clear p pref post
clear EDGELETLABPATH MATLABVERSION PATHNAMESEPARATOR
clear Friend PREFERIMAGEGRAPHICS MATLABPATHSEPARATOR

disp('Done!')


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: EdgeletPath.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
