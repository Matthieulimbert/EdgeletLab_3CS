function binary = DisplayEdgeletPath_onBinaryImage(chirpletpath,cgraphparam, sample)
% DisplayChirpletPath -- Put the pixels corresponding to the chirplet path to 1 on a binary image
%            
%  Usage
%    DisplayEdgeletPath(chirpletpath,cgraphparam, sample)
%  Inputs
%    chirpletpath   a chirplet path as returned by the function
%                   CalculateTestStatistic
%    cgraphparam    chirplet graph parameters as returned by
%                   GetChirpletGraphParam. Has to be the same set of
%                   parameters as used in CalculateTestStatistic.   
%    sample         the array containing the column indices that we sample. 
%
%  Output
%    binaryImage    a binary image whose pixel = 1 if there is an edge, and
%                   pixel = 0 if there is no edge
%
%  See Also
%    - GetChirpletGraphParam
%    - CalculateStatistic
%    - Path2TF, DisplayChirplets



%% We get the parameters :
J = cgraphparam{1}(2);    % dyadic length of signal
fsc = cgraphparam{2}(2);  % finest scale in the chirplet graph 
ymin = cgraphparam{4}(1); % lowest frequency in the chirplet graph 
ymax = cgraphparam{4}(2); % highest frequency in the chirplet graph 
delta_p = cgraphparam{5}; % Vertical distance between two ends of edgelets on the same vertical line


%% we compute the coordinates on x_axis and y-axis of :
[x,y] = Path2TF(chirpletpath,J,ymin,ymax,fsc, delta_p);    % get x y  in the sample image

x(2:end) = sample(x(2:end)); % we find the x indices in the original image

x = x + 0.5; % we change the indices to plot the edgelets on the image...
y = y + 0.5; % we change the indices to plot the edgelets on the image...


%% we plot the edgelets that are part of the path:
binary = DisplayEdgelets_onBinaryImage(x,y,2^J); % plot


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: DisplayEdgeletPath.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego

