function [costpath,edgeletpath,no_path] = CalculateStatistic_edgelet(cnetwork,maxLength, cgraphparam)
% CalculateStatistic -- for solving the optimization problems on the
%                       edgelet graphs (in C or Matlab)
%		
%  Usage
%    [costpath,edgeletpath,nedgelets] = CalculateStatistic(cnetwork,maxLength,cgraphparam)
%  Inputs
%    cnetwork       a edgelet network as returned by GetEdgeletNetwork
%    maxLength      corresponds to the maximum number of edgelets 
%                   allowed in the path. Calculates the best paths
%                   with number of edgelets equal 1,2 up to maxLength
%    cgraphparam    edgelet graph parameters as returned by GetChirpletGraphParam
%  Outputs
%    costpath       This is an array of length maxLength. Entry k in the
%                   array is the cost of the best path with k edgelets.
%    edgeletpath    This is a 1d cell array of length maxLength. Entry k in the
%                   cell array corresponds to the best path with k edgelets.
%    no_path        This is an array of length maxLength. Entry k in the
%                   array is 0 if a path with k edgelets exists, and 1 if a 
%                   path with k edgelets does not exist.
%
%  Description
%    For solving optimization longest path optimization problem
%    Before doing the optimization, the edgelet graph parameters have to set,
%    edgelet transform of data has to be done using the function and the edgelet
%    costs filled in the edgelet network.
%    Currently there is only support for edgelet graphs with connected offsets.
%    That is, near-continuity of the edgelet paths in the image
%
%  See also 
%    - EdgeletTransform for doing edgelet transforms of the data
%    - GetEdgeletNetwork for initializing the edgelet network
%    - DisplayEdgeletPath for plotting edgelet paths





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We check out whether we solve the BP optimization problem in Matlab or in C:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MEXROOT = strcat('mex',filesep,'src',filesep);
mexFiles = {... 
'BP_algo_one_path' strcat(MEXROOT,'Networks') 'cpp';
};

% check out if mex-file is installed and exists on MATLAB's search path. If not, we run the Matlab file:
if exist(mexFiles{1,1})~=3,
    type = 'M'; % we solve the optimization problem in Matlab...
else 
    type = 'C'; % We solve the optimization problem in C...
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We solve the optimization problem:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if type == 'M' % We solve the optimization problem in Matlab...
    [costpath, edgeletpath,no_path] = Solving_BP_M(...  
                                    cnetwork{1},cnetwork{2},...                                
                                    cnetwork{3},cnetwork{4},maxLength, cgraphparam);
    disp('Optimization problem resolved in Matlab');

elseif type == 'C' % We solve the optimization problem in C...
    [costpath, edgeletpath,no_path] = Solving_BP_C(...  
                                cnetwork{1},cnetwork{2},...                                
                                cnetwork{3},cnetwork{4},maxLength, cgraphparam);
    disp('Optimization problem resolved in C language');

end

% $RCSfile: CalculateStatistic.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
