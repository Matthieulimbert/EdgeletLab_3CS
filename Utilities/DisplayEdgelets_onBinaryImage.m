function binary = DisplayEdgelets_onBinaryImage(x,y,n)
% DisplayChirplets -- Put the pixels corresponding to the chirplets to 1 on a binary image
%                     
%            
%  Usage
%    DisplayChirplets(t,freq,n)
%  Inputs
%    x       array of coordinates of the edgelet path on the x-axis
%    y       array of coordinates of the edgelet path on the y-axis
%    n		 length of underlying signal
%
%  Description
%    This is a subroutine of DisplayChirpletPath so
%    it is better to use that function for plotting
%    chirplet paths. Otherwise, use Path2TF to get
%    the time and frequency coordinates for the chirplet
%    path and then use this function. 
%    
%  See Also
%    - DisplayChirpletPath, Path2TF
%    - CalculateStatistic

%% we initialize the binary image containing the path :

binary = zeros(n,n);

%% we put to 1 the pixels that correspond with the edgelet path :

nedgelets = length(x) - 1;
size_previous_edgelet = 0;

for k = 1:nedgelets,
    
  x0 = x(k); x1 = x(k+1);
  y0 = y(k); y1 = y(k+1);
  
  slope = (y1 - y0) / (x1 - x0);
  
  xx = x0 : x1 ;
  xx(length(xx)) = [];
  xx = xx - 0.5 ; % corresponds with the columns of the pixels that correspond with an edge 
  
  yy = y0 + slope * (xx - size_previous_edgelet) ; % corresponds with the lines of the pixels that correspond with an edge
  yy = yy - 0.5;
  yy = round(yy);
  
  binary ( [ yy + xx * n ]) = 1; % We put to 1 the pixels corresponding with the edgelet k
  
  size_previous_edgelet = size_previous_edgelet + length( xx ); % we update the size of the precedent edgelet
  
end

hold off;

% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: DisplayEdgelets.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego

