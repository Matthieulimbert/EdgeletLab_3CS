function cnetwork = GetEdgeletNetwork(C,cgraphparam)
% GetEdgeletNetwork -- fills in arc costs into the edgelet graph using
%                      the edgelet coefficients
%
%  Usage
%    cnetwork = GetEdgeletNetwork(C,cgraphparam)
%
%  Inputs
%    C              edgelet coefficient table as returned by EdgeletTransform.
%    cgraphparam    edgelet graph parameters as returned by GetEdgeletGraphParam. 
%                   If omitted or set to [], the default value as returned by 
%                   GetEdgeletGraphParam is used. 
%  Outputs
%    cnetwork       a cell structure of size 1-by-5, where the cell elements
%                   are G, nvertices, startnodes, endnodes, ord 
%                   in this order and are described below:
%
%        - G     Cell array of size (k,2) where k is the number of
%                nodes; cell (m,1) is a vector with all the nodes that
%                node m connects to and cell (m,2) is a vector with the 
%                corresponding arc costs
%
%        - nvertices
%                number of nodes in the network.
%
%        - startnodes  
%                List of nodes from where paths can start in the network.     
%                They must be the first nodes in the topological ordering. 
%
%        - endnodes    
%                List of nodes from where paths can end in the network.
%                They must be the last nodes in the topological ordering.
%
%        - ord	 Topological ordering of the nodes
%
%  Description
%    Preprocessing for the optimization routines in CalculateStatistic.
%    Fills in edgelet costs into the edgelet graph.
%
%  See Also
%    - EdgeletTransform
%    - GetedgeletGraphParam
%    - CalculateStatistic


% NB : this function assumes that edgelets are the arcs of the
% network. The nodes in the network are the endpoints
% of the edgelets.






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting the parameters :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We get the parameters:
n             = cgraphparam{1}(1);% size of the image
J             = cgraphparam{1}(2);% dyadic length of image
c             = cgraphparam{1}(3); % number of columns that we sample
J_c           = cgraphparam{1}(4); % dyadic length of image after we sample it
coarsestScale = cgraphparam{2}(1);  % coarsest scale in the chirplet graph 
finestScale   = cgraphparam{2}(2); % finest scale in the chirplet graph
fmin          = cgraphparam{4}(1); % lowest frequency in the chirplet graph 
fmax          = cgraphparam{4}(2); % highest frequency in the chirplet graph 
delta_p       = cgraphparam{5}; % Vertical distance between two ends of edgelets on the same vertical line


% We make sure there is no absurd value :
if (coarsestScale < 0)
    msg = strcat('The parameter coarsestScale cannot correspond',...
          ' to a negative value');
    error(msg);
end
if (finestScale < 0)
    msg = strcat('The parameter finestScale cannot correspond',...
          ' to a negative value');
    error(msg);
end
if (coarsestScale > finestScale)
    msg = strcat('The parameter coarsestScale cannot correspond',...
          ' to a value highter than the parameters finestScale');
    error(msg);
end

if (fmin < 0)
    msg = strcat('The parameter fmin cannot correspond',...
          ' to a negative value');
    error(msg);
end
if (fmax < 0)
    msg = strcat('The parameter fmax cannot correspond',...
          ' to a negative value');
    error(msg);
end
if (fmax < fmin)
    msg = strcat('The parameter fmin cannot correspond',...
          ' to a value highter than the parameters fmax');
    error(msg);
end

if (delta_p < 0 || delta_p > 1)
    msg = strcat('The parameter delta_p must of type:',...
          ' delta_p = 1 / 2^k, with k an integer');
    error(msg);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting the parameters for the best path algorithm resolution:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnetwork = cell(1,5);
nfreqs = ( fmax-fmin+1 ) / delta_p; % number of frequency indices
startfreqs = [ fmin : 1 : ( fmax-fmin+1 ) / delta_p + fmin - 1 ]; % the frequency offsets (starting frequencies)
ntfnodes = nfreqs * 2^finestScale / 2^(J-J_c); % number of left time points is 2^finestScale               
cnetwork{1} = cell(ntfnodes,2);

offset    = ( fmin /delta_p +1 : 1 : (fmax + 1) / delta_p  )';


for s = finestScale:-1:coarsestScale,
  % Calculate the changes in frequency. A time interval at scale s
  % is of length 2^(J-s)
  
  dfreq = (cgraphparam{3}{s-cgraphparam{2}(1)+1})*2^(J-s) / delta_p;
  Dfreq = repmat(dfreq,length(startfreqs),1);
  M = repmat(startfreqs',1,length(cgraphparam{3}{s-cgraphparam{2}(1)+1}));

  endfreqs = round(M+Dfreq);	% Calculate all the end frequencies and round to 
                                % get the frequency as an integer
    
  endfreqs = mod(endfreqs, n / delta_p);
  [a1 , a2] = size(endfreqs);
    
  % We fill in the vector slopeindices that contains the allowed slope indices
  for m=1:nfreqs,
      if m < (a2 + 1)/2 && m > a1 - (a2 - 1)/2
        slopeindices{m} = (a2 + 1) / 2 - m + 1 : (a2 + 1) / 2 + (a1 - m);
      elseif m < (a2 + 1)/2
        slopeindices{m} = (a2 + 1) / 2 - m + 1 : a2;
      elseif m > a1 - (a2 - 1)/2
        slopeindices{m} = 1 : (a2 + 1) / 2 + (a1 - m);
      else
        slopeindices{m} = 1:a2;
      end
  end
    
  
  % loop over time indices and find the connecting edgelets
  for b = 0:(2^s-1) / 2^(J-J_c)
    c = C{node(s,b)};
    leftoffset = b*2^(finestScale-s)*nfreqs + 1-fmin;
    rightoffset = (b+1)*2^(finestScale-s)*nfreqs + 1-fmin;
    leftindex = startfreqs + leftoffset;
    rightindex = endfreqs + rightoffset;  
    
    for m = 1:nfreqs,
       cnetwork{1}{leftindex(m),1} =... 
                  [cnetwork{1}{leftindex(m),1} rightindex(m,slopeindices{m})]; % We add the new nodes connecting the node m a each loop on each dyadic b at scale s
       cnetwork{1}{leftindex(m),2} =...
          [cnetwork{1}{leftindex(m),2} c( offset(m), slopeindices{m})];
    end
  end
end

% list the start nodes and end nodes
cnetwork{3} = 1:nfreqs; % the topological ordering of startnodes
cnetwork{4} = ntfnodes+(1:nfreqs); % the topological ordering of endnodes
cnetwork{2} = ntfnodes+nfreqs; % number of nodes in the graph

% The topological ordering (the nodes where numbered that way)
cnetwork{5} = 1:(ntfnodes+nfreqs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We write the data in files in order to test the BP algorithm resolution in C language :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for k = 1 : length(cnetwork{1})
%         destination = 'C:\Users\Matthieu\Desktop\Manip_project2\EdgeletLab_3C\Data_C\';
%         name_nodes = strcat( 'nodes_' , num2str(k) );
%         name_costs = strcat( 'costs_' , num2str(k) );
%         
%         fileID_nodes = fopen(strcat(destination,name_nodes,'.txt'),'wt');
%         fileID_costs = fopen(strcat(destination,name_costs,'.txt'),'wt');
%         
%         nodes = ( cnetwork{1}{k,1} )';
%         costs = ( cnetwork{1}{k,2} )';
%                
%         for i = 1 : length(nodes)
%             fprintf(fileID_nodes , num2str(nodes(i)) ); % we writes the nodes that goes out from node k;
%             fprintf(fileID_costs , num2str(costs(i)) ); % we writes the costs of the nodes that goes out from node k;
%             fprintf(fileID_nodes,'\n');
%             fprintf(fileID_costs,'\n');
%         end
%         fclose(fileID_nodes);
%         fclose(fileID_costs);
% end

% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: GetEdgeletNetwork.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego

