function C = ET_standard(I_b, param)

% ET_standard -- Calculates edgelet coefficients 
%            
%  Usage
%    C = ET_standard(I_b, param, delta_p);
%  Inputs
%    I_b     Noisy image of size n x n pixels
%    param   edgelet graph parameters as returned by GetChirpletGraphParam
%  Outputs
%    C       edgelet coefficients, stored in a cell array.
%
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
%    - EdgeletTransform
%    - GetEdgeletGraphParam for information on how to customize the chirplet transform.
%    - GetEdgeletNetwork, CalculateStatistic, for setting up the 
%      edgelet graph and calculating the test statistic. 




%% Getting some parameters :

n   = param{1}(1);
J   = param{1}(2);
c   = param{1}(3);
J_c = param{1}(4);
delta_p = param{5}; % Vertical distance between two ends of edgelets on the same vertical line
coarsestScale = param{2}(1);
finestScale = param{2}(2);
w = param{7}; % Depth of the window filter:

nCells = 2^(finestScale+1) - 1; % would need to change node.m if the stuff in previous line is used
C = cell(1,nCells);

% Minimum and maximum values of the offsets (in order to not compute the filter close to the borders on the image (on the top and on the bottom) :
minfreq = param{4}(1);
maxfreq = param{4}(2);

% vector that contains the values of the offsets :
offset    = ( minfreq /delta_p +1 : 1 : (maxfreq + 1) / delta_p - 1 )';
nb_offset = length(offset);


%% We transform I_b into a sub-pixelique image :
I_b_sub       = repmat (I_b, 1 ,1 / delta_p);
I_b_sub       = reshape(I_b_sub', c, n / delta_p)'; % On ne fait que redonner la même valeur aux sous-pixels ici... On pourrait faire une interpolation sous pixelique...

[size_x , size_y ] = size(I_b_sub);


%% Convolution of all of the image :
I_b_conv  = I_b_sub;
I_b_conv  = I_b_conv(:);
vect_conv = - [ ones(w / delta_p , 1); - ones(w / delta_p , 1)];
I_b_conv  = conv(I_b_conv, vect_conv, 'same');
I_b_conv  = reshape(I_b_conv, n / delta_p, c );

I_b_conv(1: w / delta_p , : )                   = 0; % we don't want to include the edges of the image in the convolution
I_b_conv( (n - w) / delta_p : n / delta_p , : ) = 0; 





%% Loop For on the scale : 2^-finestScale to 2^-coarsestScale :
for s = finestScale:-1:coarsestScale
    
  slopes = param{3}{s-coarsestScale+1};
  nSlopes = length(slopes);
  

    
    %% Loop For on the dyadic at the scale 2^-s :
    for b = 0:(2^s-1) / 2^(J-J_c)
        

        % initialization of the matrix containing the coefficient for the
        % dyadic b at the scale s :
        mat_rep_filter = zeros( n / delta_p , nSlopes);
        
        % handle dyadic interval (s,b)
        ix = dyadindex(s,b,n);
        size_dyadic = length (ix);     
        
        % Convolution of the image in the dyadic considered :
        I_b_conv_temp = I_b_conv ( : , ix(1) : ix(size_dyadic) );
        pos = [0.5 : 1 : size_dyadic - 0.5];
        v_1 = 0 : size_dyadic - 1;
        
        
        %% First methods with loops "for" on the offsets and on the slopes :
        for i = 1 : nb_offset;
            for j = 1 : nSlopes
                
            pos_edgelet = offset(i) * delta_p + pos * slopes(j);
            pos_edgelet = round(pos_edgelet) ;
            v = pos_edgelet / delta_p + size_x * v_1;
            
            if max(v) > size_dyadic * size_x
               rep_filter = 0; 
            else
               rep_filter = I_b_conv_temp ( v );
               % rep_filter = sum( rep_filter ) / (size_dyadic * w / delta_p);
               rep_filter = sum( rep_filter ) / sqrt(2 * (size_dyadic * w / delta_p) );
            end
            
            mat_rep_filter (offset(i),j) = rep_filter;
            
            
            end 
        end
        
        
         %% Second method without the loops "For" on the offsets and the slopes (greedy in termps of memory...) :
         
%         offset_1  = repmat(offset, 1, size_dyadic * nSlopes);
%         pos_1     = repmat(pos,  length(offset), nSlopes );
%         slopes_1  = slopes(:);
%         slopes_1  = repmat(slopes_1, 1, size_dyadic);
%         slopes_1  = reshape(slopes_1', 1, size_dyadic * nSlopes);
%         slopes_1  = repmat(slopes_1, length(offset), 1);
%         
%         pos_edgelet = offset_1 * delta_p + pos_1 .* slopes_1;
%         pos_edgelet = round(pos_edgelet) ;
%         v_1 = repmat(v_1, length(offset), nSlopes);
%         v = pos_edgelet / delta_p + size_x * v_1;
%         
%         v( v > size_dyadic * size_x ) = 1; % we remove the absurd values
%         rep_filter = I_b_conv_temp ( v(:) );
%         rep_filter = reshape(rep_filter, length(offset), size_dyadic * nSlopes);
%         
%         rep_filter = rep_filter';
%         rep_filter = reshape(rep_filter(:), size_dyadic , nSlopes * length(offset) );
%         rep_filter = sum(rep_filter) / sqrt(2 * (size_dyadic * w / delta_p) );
%         rep_filter = reshape(rep_filter, nSlopes, length(offset) )';
%         
%         mat_rep_filter_2 = [ zeros( offset(1) - 1, nSlopes) ; rep_filter ; zeros(n / delta_p - (offset(1) - 1 + length(offset)) , nSlopes) ];
        
        
        

        
        %% We put in memory the coefficients of the edgelets in the dyadic b at the scale s :
        C{node(s,b)} = mat_rep_filter ; % Each edgelet has a weight in relation with its length...

        
    end
    
end


% $RCSfile: ET_standard.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego
