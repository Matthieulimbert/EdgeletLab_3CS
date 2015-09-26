function sample = sampling( cgraphparam )
% sample --  Returns the x indices of the columns of the image that we
%            sample. If c = n, we take all the columns of the image
%  Usage
%   [t,freq] = Path2TF(p,J,fmin,fmax,fsc, delta_p)
%  Inputs
%    cgraphparam		 edgelet graph parameters as returned by GetEdgeletGraphParam. 
%                        If omitted or set to [], the default value as returned by 
%                        GetEdgeletGraphParam is used.    
%  Outputs
%    sample              array containing the x indices of the columns that
%                        we sample
%
%  See Also
%    - DisplayEdgelets, DisplayEdgeletPath
%    - CalculateStatistic
%    - Solving_BP_M, Solving_BP_C





%% Getting the parameters :
n   = cgraphparam{1}(1);
J   = cgraphparam{1}(2);
c   = cgraphparam{1}(3);
J_c = cgraphparam{1}(4);



%% Getting the column indices :
if J == J_c
   sample = 1:n ; % we consider all the pixels 
else
    r = (n-2) / (c-2);
    vect = 1 : c-2;
    vect = floor(vect * r);
    sample = horzcat( 1 , vect , n); %we consider only c columns of pixels
end


end


% $RCSfile: sampling.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
