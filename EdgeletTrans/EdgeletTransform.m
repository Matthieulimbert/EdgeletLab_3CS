function C = EdgeletTransform(I_b , param)
% EdgeletTransform
%
%  Usage
%    C = EdgeletTransform(I_b,param)
%  Inputs
%    I_b     Noisy image of size n x n pixels
%    param   edgelet graph parameters as returned by GetEdgeletGraphParam
%  Outputs
%    C       edgelet coefficients, stored in a 1d cell array. 
%            The length of the array is equal to the number of dyadic time intervals.
%
%  Description
%    This function computes the edgelet responses of all the edgelets in
%    the edgelet dictionnary.
%
%    The data structure used for storing the edgelet coefficients is a Matlab cell array.
%    Each element in the array is a table of edgelet coefficients corresponding 
%    to a dyadic time interval [k 2^(-j), (k+1)2^(-j)], labeled (j,k),
%    where k=0,1,...,2^j-1 and j is a scale index which can range from 0 to J-1.
%    A table corresponding to a dyadic interval (j,k) is indexed by (2^j + k)
%    and is an n-by-(#slopes at scale j) matrix whose values are described as follows:
%      -each column in the table corresponds to one slope,
%       column k corresponds to the k-th slope
%      -each entry in a column corresponds to a frequency offset,
%       the m-th entry in each column corresponds to the frequency offset 2*pi*(m-1)/n
%
%  See Also
%    - GetEdgeletGraphParam for information on how to customize the chirplet transform.
%    - GetEdgeletNetwork, CalculateStatistic, for setting up the 
%      edgelet graph and calculating the test statistic. 

	C = ET_standard(I_b,param);
    
end

% $RCSfile: EdgeletTransform.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego


