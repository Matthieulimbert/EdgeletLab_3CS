function slopes = GetSlopes_edgelet(parameters,sc)
% GetSlopes -- Returns a list of chirplet slopes
%
%  Usage
%    slopes = GetSlopes(parameters,sc);
%    slopes = GetSlopes(parameters);
%  Inputs
%    parameters  a cell array with parameters required for the slope
%                discretization scheme. Entries should be as follows:
%      parameters{1}  vector of length 2 or 
%                     a matrix of size 2-by-(num. scales).
%                     - In the first case, sets the range of slopes to 
%                       [parameters{2}(1),parameters{2}(2)].
%                     - In the latter case, sets the range of slopes for
%                       each individual scale, from the coarsest to
%                       the finest scale. The first row is for the
%                       lower bound on the slope, the second row for
%                       the higher bound.
%      parameters{2}  vector with discretization step of slopes at a 
%                     particular scale, going from the coarsest to
%                     the finest scale.
%      parameters{3}  vector of length 2, first entry is the coarsest 
%                     allowable scale, second entry is the finest
%                     allowable scale.
%    sc          to be used to return only the slopes at the particular 
%                scale 2^(-sc).
%  Outputs
%    slopes      The slopes set by the parameters.
%                If sc was set, it is a an array, 
%                otherwise it is a cell array of size 1-by-(#scales), with 
%                entry k corresponding to the k-th coarsest scale. 
%
%
%  See Also:
%     - GetChirpletGraphParam

if (nargin<2),
    nScales = parameters{3}(2)-parameters{3}(1)+1;
    slopes = cell(1,nScales);
    for k=1:nScales,
        if (isvector(parameters{1})),
            % same slope range for every scale
            slopes{k} = (parameters{1}(1):parameters{2}(k):parameters{1}(2)); 
        else
            % different slope range for every scale
            slopes{k} =... 
                (parameters{1}(1,k):parameters{2}(k):parameters{1}(2,k)); 
        end
    end
else
    ind = parameters{3}(1)+1-sc;
    slopes = parameters{1}(1):parameters{2}(ind):parameters{1}(2); 
end

% $RCSfile: GetSlopes.m,v $
% $Date: 2006/05/01 17:44:15 $
% $Revision: 1.3 $
%
% Copyright (c) Hannes Helgason, California Institute of Technology, 2006
