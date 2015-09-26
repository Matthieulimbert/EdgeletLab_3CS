function CompileMex
% CompileMex -- checks whether mex-files have been compiled and compiles them if needed

MEXCOMPILED = 1;

MEXROOT = strcat('mex',filesep,'src',filesep);

% The list of mex-files
% Every row in the string array should be
%'NAME-OF-FILE' 'LOCATION-RELATIVE-TO-CHIRPLABROOT' 'FILE-EXTENSION'


mexFiles = {... 
'BP_algo_one_path' strcat(MEXROOT,'Networks') 'cpp';
'csa' strcat(MEXROOT,'Networks') 'cpp';
};

nFiles = size(mexFiles,1);

% Check if MEX files are installed
for k=1:nFiles,
  % check if mex-file exists on MATLAB's search path
  if exist(mexFiles{k,1})~=3,
    MEXCOMPILED = 0;
    break;
  end
end

% If MEX files are not installed, try to install them
if ~MEXCOMPILED,
  disp('MEX files are not installed')
  disp('Compiling MEX files...')
  
  rootDir = pwd;
  for k=1:nFiles,
    cd(mexFiles{k,2});
    fullFileName = sprintf('%s%s%s.%s',mexFiles{k,2},filesep,mexFiles{k,1},mexFiles{k,3});
    fprintf('\n'); % go to line
    fprintf('Compiling: %s',fullFileName);
    fprintf('\n'); % go to line
    try,
      eval(sprintf('mex %s.%s',mexFiles{k,1},mexFiles{k,3}));
    catch
      fprintf('CompileMex: Error occured while compiling MEX file: %s',fullFileName)
      disp('Try running "mex -help" to see if mex compiler is properly installed.')
      disp('If problem persists, EdgeletLab will try to run a slower m-file versions of this function. ')
      % go back to EdgeletLab root directory
      cd(rootDir);
    end
    % go back to EdgeletLab root directory
    cd(rootDir);
  end
  fprintf('\n'); % go to line
else
     disp('MEX files are already installed')
end


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: CompileMex.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
