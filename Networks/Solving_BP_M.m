function [cost,paths,no_path] = Solving_BP_M(G,n,startnodes,endnodes,maxLength,cgraphparam)        
% Solving_BP_M --   solves the longest path optuimization problem for an acyclic graph with
%                   a constraint on the maximum number of arcs that the path can have. This
%                   script solves the problem entirely in Matlab.
%
%  Usage
%    [cost,paths,no_path] = Solving_BP_M(G,n,startnodes,endnodes,maxLength,cgraphparam)
%  Inputs
%    G              a cell array of size (k,2) where k is the number of nodes in the graph
%                   for which arcs go out from. cell (m,1) is a vector with all the nodes 
%                   that node m connects to and cell (m,2) is a vector with the corresponding
%                   arc costs
%    n              total number of nodes in graph.
%    startnodes     list of nodes from where paths can start in the acyclic network.
%    endnodes       list of nodes from where paths can end in the acyclic network.
%    maxLength      maximum number of arcs that the best path is allowed to have.
%    cgraphparam    edgelet graph parameters as returned by GetChirpletGraphParam
%
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
%    This is an inner function in EdgleletLab.
%    It is recommended to use it through CalculateStatistic.
%
%  See Also:
%    - CalculateStatistic
%    - Solving_BP_C




%% We get the parameters :

J       = cgraphparam{1}(2);    % dyadic length of signal
fsc     = cgraphparam{2}(2);    % finest scale in the chirplet graph 
fmin    = cgraphparam{4}(1);    % lowest frequency in the chirplet graph 
fmax    = cgraphparam{4}(2);    % highest frequency in the chirplet graph 
delta_p = cgraphparam{5};       % Vertical distance between two ends of edgelets on the same vertical line (in pixels)
C_alpha = cgraphparam{6}(1);    % parameter of the horizon curve
alpha   = cgraphparam{6}(2);    % parameter of the horizon curve

angle_max = pi/8; % depends on C_alpha, alpha, and the scale...

nG = size(G,1);  % returns the number of starting nodes in G
% Assume that the nodes in G appear in topological order.
ord = 1:nG;



%% We initialize the matrices d and pred:

% Initialize the distance label matrix   
d = zeros(maxLength+1,n); % we want to maximize the path, so all the startnodes have distance label 0 

% Initialize the predecessor indices
pred = zeros(maxLength+1,n);     % we initialize the predecessor labels to zero
pred(1,startnodes) = startnodes; % we create predecessors for startnodes and l = 1



%% Examine nodes in topological order and get the distances and predecessors of each node:
for start = 1:nG,
    A = G{start,1}; % A is the set of endpoints and costs of arcs emanating out of node 'start'.
    L = G{start,2}; % costs

    for k = 1:length(A), % we look at all the arcs emanating out of node 'start'.

        for l=2:(maxLength+1), % we look at those arcs for all the lengths allowed (l= 1 : maxlength). 

            if pred(l-1,start) == 0 % we check out if the path that we are computing makes sense..
                continue;
            end

            % We check out whether the angle between the node start and the node A(k) is too large or not:
            if start > length(startnodes)
                [t,f] = Path2TF( [ pred(l-1,start), start , A(k)] , J , fmin , fmax , fsc , delta_p );    % get time-freq indices
                a = ( f(2) - f(1) ) / ( t(2) - t(1) );
                b = ( f(3) - f(2) ) / ( t(3) - t(2) );
                angle = atan( (a-b) / (1 + a*b) );
                if angle > angle_max
                    continue;
                end
            end

            if ( abs( d(l,A(k)) ) <= abs( ( d(l-1,start)+L(k) ) ) )
                d(l,A(k)) = d(l-1,start)+L(k);
                pred(l,A(k)) = start;
            end
        end
    end
end

  

%% We find the costs and the endnodes with highest distances :

% Now we figure out which of the endnodes has the hightest distance label. That is going to be
% the end node of our longest path. Then, to find the longestest path we just have to start at
% the endnode and track the predecessor nodes.
cost = zeros(1,maxLength+1);
ind = zeros(1,maxLength+1);

for m = 1:(maxLength+1),
    [cost(m),ind(m)] = max( abs (d(m,endnodes)) ); % we maximize the cost divided by the length of the path, in absolut value
end
cost = cost(2:(maxLength+1));



%% We fin the paths of sizes l = 1, ..., maxLength:
paths = cell(1,maxLength);
no_path = zeros(1,maxLength); % vector that indicates if there is a path or not for each length (0 means there is a path, 1 means that there is no path)

for k=1:maxLength,
paths{k} = zeros(1,k+1);
paths{k}(k+1) = endnodes(ind(k+1));

for l=k:-1:1,
  if paths{k}(l+1) == 0 % There is not any path...
    paths{k}(l) = 0; 
    no_path(k) = 1;
  else
    paths{k}(l) = pred(l+1,paths{k}(l+1));
    if paths{k}(l) == 0 % There is not any path...
        no_path(k) = 1;
    end        
  end
end

end
  
  
end




% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: Main_EdgeletLab.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
