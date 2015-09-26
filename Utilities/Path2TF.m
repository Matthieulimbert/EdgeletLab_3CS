function [x,y] = Path2TF(p,J,ymin,ymax,fsc, delta_p)
% Path2TF -- Returns the x and y indices
%            for the given edgelet path p as returned by
%            CalculateStatistic, minctratiocells or ShortestPathCell.
%  Usage
%   [t,freq] = Path2TF(p,J,fmin,fmax,fsc, delta_p)
%  Inputs
%    p		 array of the topological order of the nodes along a path
%    J       dyadic length of signal
%    ymin    minimum y coordinate from which pixels are computed in the image in EdgeletLab.
%            (warning : y-axis starts on the upper edge of the image and points
%            to the bottom on the image...), default is 0
%    ymax    maximum y coordinate until which pixels are computed in the image in EdgeletLab. 
%            (warning : y-axis starts on the upper edge of the image and points 
%            to the bottom of the image...), default is 2^J - 1
%    fsc     the finest scale in the chirplet gragh, default J-1
%  Outputs
%    x       array of coordinates of the edgelet path on the x-axis
%    y       array of coordinates of the edgelet path on the y-axis
%
%  See Also
%    - DisplayEdgelets, DisplayEdgeletPath
%    - CalculateStatistic
%    - Solving_BP_M, Solving_BP_C


%% We check the inputs:

if (nargin < 6),
    delta_p = 1;
    if (nargin < 5)
        fsc = J-1;
        if (nargin < 4),
            ymax = 2^J-1;
            if (nargin < 3),
                ymin = 0;
            end
        end
    end
end


%% We compute the coordinates of the path:

nfreqs = ( ymax-ymin+1 ) / delta_p; % number of y indices per x
 
y = mod(p-1,nfreqs);

x = (p - y - 1)/(nfreqs/2^(J-fsc)); % PEND! Needs to be changed for finest scale         

y = y * delta_p; % We need to change the value of the y coordinate if we are considering an sub-pixelic image

y = y + ymin + 1 * delta_p;

             

             

% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: Path2TF.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
