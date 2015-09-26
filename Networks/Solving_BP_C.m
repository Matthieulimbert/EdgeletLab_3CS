function [cost,paths, no_path] = Solving_BP_C(G,n,startnodes,endnodes,maxLength, cgraphparam) 

% Solving_BP_C --   solves the longest path optimization problem for an acyclic graph with
%                   a constraint on the maximum number of arcs that the path can have. This
%                   script solves the problem calling a mex function in C language.
%
%  Usage
%    [cost,paths,no_path] = Solving_BP_M(G,n,startnodes,endnodes,maxLength,cgraphparam)
%
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
%    - Solving_BP_M




%% Call mex function to examine nodes in topological order and get the distances and predecessors of each node:
[d,pred] = BP_algo_one_path(G(:,1),G(:,2),n,startnodes,endnodes, maxLength, endnodes(1), cgraphparam);



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
if (nargout > 1),
  % the caller wanted to get the paths, so find them 
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


% $RCSfile: Solving_BP_C.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego

