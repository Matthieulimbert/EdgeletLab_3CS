function DisplayEdgelets(x,y,n)
% DisplayChirplets -- Plots an edgelet path given 
%                     x-axis and y-axis indices
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



nedgelets = length(x) - 1;

for k = 1:nedgelets,
    
  x0 = x(k); x1 = x(k+1);
  y0 = y(k); y1 = y(k+1);
  
  slope = (y1 - y0) / (x1 - x0);
  xx = [x0 x1];
  
  omega = y0 +  slope .* (xx - x0); 
  
  plot([x0 x0], [0.5 n+0.5],'r--'); plot([x1 x1], [0.5 n+0.5],'r--'); % We change thse scale fot plotting on the image
  plot(xx,omega,'LineWidth',2,'Color','r');
  
end

hold off;

% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: DisplayEdgelets.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego

