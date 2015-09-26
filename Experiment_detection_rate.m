%
% Main script - EdgeletLab_3CS
% You must run EdgeletPath.m before.
%
%
% Steps:
% 1) Initialization of the parameters (image and graph)
% 2) Generate an image
% 3) Sample the columns of the image (if needed)
% 4) Calculate the edgelet responses (Edgelet Transform)
% 5) Fill in edgelet costs in network
% 6) Calculate the Best Path for prescribed path lengths test statistic
% 7) Plot edgelet paths corresponding to the lengths chosen for the BP test statistic.
%
% You must run EdgeletPath.m before.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cleaning the workplace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clear path

disp('--Running Main_EdgeletLab.m--');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION OF THE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

  %WE DO NOT RUN AN EXISTING SCENARIO AND GIVE THE PARAMETERS MANUALLY:

    % We manually enter the parameters through the Main_EdgeletLab.m...   

    % Size of the image :
    J = 7;                                      % dyadic size of image
    n =  2^J;                                   % size of image (Note that the image is a square).
    
    % Sampling parameters :
    J_c = 7;                                    %(must be <= J)
    c   = 2^J_c;                                % number of columns we sample

    % Mean and variance of the additive gaussian noise:
    M = 0;                                      % mean
    %V = 1;                                      % variance :

    % Vertical distance between two ends of edgelets on the same vertical line :
    delta_p = 1;                                % (in pixel) It must be of the form : delta_p = 1 / 2^k, with k is an integer
    delta = delta_p / n;                        % (in distance)

    % Maximum slope of an edgelet (in absolut value) :
    slopeRange = 1;                             % (45 degrees slope, if we consider the entire image)
    slopeRange = slopeRange * 2^(J - J_c);      % we adjust the slope considering the sampling 

    % depth of the filter :
    w = 3;                                      % (in pixels)

    % Upper bound and lower bound of the horizon function on the y-axis (in pixels) to not compute all the image :
    ymin = w;
    ymax = n - w;


    % Coarsest and finest scales :
%     csc = J - J_c;                              % Coarsest scale we can compute (warning :must be >= J - J_c) 
% %     csc = 2;
%     fsc = J-1;                                  % finest scale we can compute 
    csc = J - J_c+3;                              % Coarsest scale we can compute (warning :must be >= J - J_c) 
    fsc = csc;                                  % finest scale we can compute 



    % Type of image :
%     type_image = 'line';
    type_image = 'sin';
%     type_image = 'noise_only';
%     type_image = 'sin_am';


    % Maximum size of the path :
    maxLength = 8;                              % Calculates BP for lengths 1,...,maxLength
    
    
    % Parameters of the horizon curve :
    C_alpha = 2;
    alpha   = 1.7;




    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE THE NOISELESS IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating an image...');
V = 0;
[I, I_b] = MakeImage(n, M, V, type_image);% I is the noiseless image, I_b is the noisy image




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCALIZING THE HORIZON FUNCTION FOR COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We perform a Canny filter on the noiseless image to get a matrix showing where the horizon function really is:

thresh            = [0.2 , 0.4];                % Canny filter threshold
I_comparison      = edge(I,'canny',thresh);   % We perform the Canny filter
I_comparison(:,1) = 1;
I_comparison(:,n) = 1;


figure;
imshow(I_comparison);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STARTING OF THE LOOP ON THE VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
for V = 0 : 2
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE THE NOISY IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating an image...');
[I, I_b] = MakeImage(n, M, V, type_image);% I is the noiseless image, I_b is the noisy image




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY THE IMAGES GENERATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Operations to be able to plot the images:
I_b_temp = I_b;
I_b_temp (I_b_temp > 1) = 1; %pixels have to be between 0 and 1.
I_b_temp (I_b_temp < 0) = 0;

I_temp = 1 - I; % In Matlab, if pixel=0, then it's black...
I_b_temp = 1 - I_b_temp;
 
% % Plot the noiseless image:
% figure(1)
% imshow(I_temp);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Taking the parameters...');
    graphparam = GetEdgeletGraphParam(n, csc, fsc, delta_p, slopeRange, ymin, ymax, C_alpha, alpha, w, c);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling the c columns of the images (if needed, i.e c < n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = sampling( graphparam );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET EDGELET TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Taking edgelet transform...');
cc = EdgeletTransform(I_b(:,sample), graphparam); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTATION OF THE BEST PATH:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preprocessing for the computation of the best path :
disp('Generating chirplet graph and assigning costs... (this will take a while)');
cnetwork = GetEdgeletNetwork(cc,graphparam);

% Computation of the best path :
disp('Running optimization routine for graph...');
[costpath, bestpaths, no_path] = CalculateStatistic_edgelet(cnetwork,maxLength, graphparam);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Getting the threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_s = 0;


% We reject the paths with responses lower than the threshold :
no_path (costpath < t_s) = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT THE BEST PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp('We display the images with the paths of each length...')

npaths = length(bestpaths);
binary_image = cell(1,npaths);

for k=1:npaths,
  if no_path(k) == 0 % If a path exists for this length, we plot it...
       h = figure(k + 2);

      % subplot(npaths,1,k);
      imshow(I_b_temp);
      hold on;
    
      binary_image{k} = DisplayEdgeletPath_onBinaryImage(bestpaths{k},graphparam, sample);
      figure();
      imshow(binary_image{k});
      
      DisplayEdgeletPath(bestpaths{k},graphparam, sample); % plot the best path

      title(strcat('Best path restricted to #edgelets=',num2str(k),...
                       ', |cost|=',num2str(abs(costpath(k)))));
                   
        

      
  else
      fprintf( 'Warning : There is no path of length l = %d \n', k);
  end
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WE CHECK IF THE EDGELET PATH DOES CORRESPOND TO THE HORIZON FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

DisplayEdgeletPath(bestpaths{k},graphparam, sample); % plot the best path



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF THE LOOP ON THE VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
end


disp('Done!');


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: Main_EdgeletLab.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego