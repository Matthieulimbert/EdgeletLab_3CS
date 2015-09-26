function parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq,maxfreq, C_alpha, alpha, w, c)
% GetEdgeletGraphParam -- parameter configuration for the edgelet graph
%
%  Usage
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq,maxfreq, C_alpha, alpha, w, c)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq,maxfreq, C_alpha, alpha, w)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq,maxfreq, C_alpha, alpha)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq,maxfreq, C_alpha)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq,maxfreq)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange,minfreq)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p,slopeRange)
%    parameters = GetEdgeletGraphParam(n,csc,fsc,delta_p)
%    parameters = GetEdgeletGraphParam(n,csc,fsc)
%    parameters = GetEdgeletGraphParam(n,csc)
%    parameters = GetEdgeletGraphParam(n)
%
%  Inputs
%    n          size of the image, n=2^J
%    csc        the coarsest scale index for which to calculate the chirplet transform,
%               corresponds to the scale 2^(-csc). 
%               Has to be a value smaller than or equal to the input parameter fsc.
%               If value is set to [], default value will be used.
%               Default value is csc=0.
%    fsc        the finest scale for which to calculate the chirplet transform,
%               corresponds to the scale 2^(-fsc).
%               Has to be a value bigger than or equal to csc.
%               If value is set to [], default value will be used.
%               Default value is fsc=J-1.
%    delta      Difference on the verticale line between the rigth
%               extremities of two edgelets of the same dyadic.
%    slopeRange Sets the range of slopes. 
%               Can either be a scalar or an array of length 2. 
%               In case of a scalar, the slope range will be set to [-slopeRange,slopeRange].
%               In case of an array, the slope range will be set to [slopeRange(1),slopeRange(2)] 
%               If omitted or value is set to [], the default range [-1,1] will be used.
%    minfreq    minimum y coordinate from which pixels are computed in the image in EdgeletLab.
%               (warning : y-axis starts on the upper edge of the image and points
%               to the bottom on the image...), default is 0
%    maxfreq    maximum y coordinate until which pixels are computed in the image in EdgeletLab. 
%               (warning : y-axis starts on the upper edge of the image and points 
%               to the bottom of the image...), default is 2^J - 1. We must have minfreq <= maxfreq <= n-1.
%               If value is set to [], default value will be used.
%               Default value is n-1.
%    C_alpha    parameter of the horizon curve (must be positive)
%    alpha      parameter of the horizon curve (we must have 1 < alpha <= 2
%    w          depth of the window filter (must be an integer > 0)
%    c          number of columns that we want to sample (if c = n, there
%               is no sampling). Must be of the form c = 2^J_c, with 0 <=
%               J_c <= J, and integer
%    
%
%  Output
%    parameters  a cell array with parameters for the chirplet graph and transform.
%                The entries are:
%      - parameters{1}  vector of length 2, first entry lenght of signal, N,
%                       second entry dyadic length of signal, J (N=2^J)
%      - parameters{2}  vector of length 2, first entry is the coarsest 
%                       scale for which to calculate the chirplet transform, 
%                       second entry is the finest scale for which to
%                       calculate the chirplet transform.
%      - parameters{3}  slopes stored in a cell array of size 
%                       1-by-(number of allowable scales).
%                       Entry k has the slopes at the k-th coarsest scale,
%                       i.e. goes from coarsest to finest scale.
%      - parameters{4}  vector of length 2, first entry is the lowest frequency
%                       in the graph, second entry is the highest frequency.
%                       1-by-(number of allowable scales).
%                       Entry k has the slopes at the k-th coarsest scale,
%                       i.e. goes from coarsest to finest scale.
%      - parameters{5}  scalar: Vertical distance between two ends of
%                       edgelets on the same vertical line. Must be of the form : 1 / 2^k,
%                       with k an integer >= 0.
%      - parameters{6}  vector of length 2 containing the parameters of the
%                       horizon curve. First entry is the parameter C_alpha, second entry is
%                       the paramerer alpha.
%      - parameters{7}  scalar containing w the depth of the window filter
%
%
%  Description
%    Creates a cell array with parameter configuration which determine the edgelet graph. 
%    These parameters are then later on passed to the edgelet transform, optimization
%    routines and plotting utilities. The main purpose of wrapping these parameters in one 
%    structure is to simplify parameter passing between functions. 
%
%    NOTE: This function only allows for evenly spaced slopes at each fixed scale but
%    it is possible to overwrite the slope discretization parameter. See the user's manual.
%
%  See Also
%    - GetSlopes - can be used to custom make slope discretization schemes

% PEND! Perhaps this structure should also have information about
%       number of edgelets in the graph, whether it assumes connected
%       offsets or slopes (would need to include the curvature constraint)
%       and so forth.
% PEND! Switch to struct if it does not affect performance. 
%       The output parameters, parameters{4} etc. look cryptic.
% PEND! For input parameters, switch to ('Property','Value') style.






%% Get the size of the image:

J   = ceil(log2(n)); % dyadic length of signal
J_c = ceil(log2(c)); % dyadic length of signal

if (2^J ~= n),
  msg = 'GetEdgeletGraphParam: image size is not dyadic. Should be of the form n=2^J.';
  msgid = 'ChirpLab:ChirpletGraph:NonDyadicLength';
  error(msgid,msg);
end


%% Check number of arguments and set to default values:

if (nargin < 11)
    c = 2^J; % default sampling : we take all the pixels of the image
    if (nargin < 10)
        w = 3;
        if (nargin < 9),
            alpha = 1.7;
            if (nargin < 8),
                C_alpha = 2;
                if (nargin < 7),
                    maxfreq = n - w;
                    if (nargin < 6),
                        minfreq = w;
                        if (nargin < 5),
                            % use default slope discretization range
                            slopeRange = 1;
                            if (nargin < 4),
                                delta_p = 1;
                                if (nargin < 3),
                                    fsc = J-1;
                                    if (nargin < 2),
                                        csc = 0;
                                    end
                                end
                            end
                        end    
                    end
                end
            end
        end
    end
end


%% Check if any parameter has an empty value
if(isempty(csc)),
  csc = 0; 
end
if(isempty(fsc)),
  fsc = J-1;
end


%% Check the value of slopeRange:
if(length(slopeRange) == 0),
  % empty slopeRange, use default
  slopeRange = 1;
elseif(length(slopeRange) == 1),
  % slopeRange is a scalar
  slopeRange = [-slopeRange slopeRange];
elseif (length(slopeRange) > 2)
  errormsg = 'GetEdgeletGraphParam: slopeRange can only be empty or scalar or vector of length 2';
  error(errormsg);
end;


%% Now check the final value of slopeRange:
if (slopeRange(1) > slopeRange(2))
  errormsg = 'GetChirpletGraphParam: slopeRange(1) must be less than slopeRange(2)';
  error(errormsg);
end;
if(isempty(minfreq)),
  minfreq = w; 
end
if(isempty(maxfreq)),
  maxfreq = n-w; 
end

if (csc>fsc),
  msg = 'ERROR: The parameters csc and fsc must satisfy csc<=fsc.';
  msgid = 'EdgeletLab:GetEdgeletGraph:ScalesNotRight';
  error(msgid,msg);
end

if (csc< J - J_c),
  csc = J - J_c; % we change the coarsest scale if it does not respect the condition : csc >= J - J_c 
end

if (minfreq > maxfreq || minfreq < 0 || maxfreq >= n),
  msg = 'ERROR: The parameters minfreq and maxfreq must satisfy 0<= minfreq <= maxfreq <= n-1.';
  msgid = 'EdgeletLab:GetEdgeletGraph:ForbiddenFrequencyRange';
  error(msgid,msg);
end


%% Set the edgelet graph parameters
parameters = cell(1,7);
parameters{1} = [n J c J_c];
parameters{2} = [csc fsc];

%% Get the slopes
slopeparam = cell(1,3);
slopeparam{1} = slopeRange;
slopeparam{2} = zeros(1,fsc-csc+1); % #allowable scales=fsc-csc+1
slopeDifference = abs(slopeRange(2)-slopeRange(1));
for s=csc:fsc,
    nx = 2^(J-s);
    sldf = delta_p * 2 / slopeDifference ;
    
%    sldf       for adjusting how fine/coarse the slope discretization is. 
%               Sets the difference between adjacent slopes at scale 2^(-s) to
%                     sldf * (maximumslope-minimumslope)/2 * 2^(s-J)
%               The smaller the value is, the finer the slope discretization will be.
%               If it is a power of 2 the time-frequency endpoints of the chirplets 
%               coincide with the discrete frequencies.
%               If value is set to [], default value will be used.
%               Default value is 2.

    if (sldf/nx > 1),
      % just to make sure that we have at least 2 slopes per scale
	  % in the case of the finest allowable scale (s=J-1)
      slopeparam{2}(s-csc+1) = slopeDifference/2;
    else
      slopeparam{2}(s-csc+1) = (slopeDifference/2)/nx*sldf; % Correspond to the definition of sldf (see above).
    end
end
slopeparam{3} = parameters{2};

parameters{3} = GetSlopes_edgelet(slopeparam);


%% set maximum and minimum frequencies
parameters{4} = [minfreq maxfreq];


%% Vertical distance between two ends of edgelets on the same vertical line :
parameters{5} = delta_p; % in pixels


%% parameters of the horizon curve:
parameters{6} = [C_alpha , alpha]; % so that we can compute the maximum angle in the BP algorithm (depending on the scale)


%% Depth of the window filter:
parameters{7} = w;


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: GetEdgeletGraphParam.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego


