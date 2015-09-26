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

% Use_scenario:
%
% - If Use_scenario = 0 : we do not use an existing scenario, and we enter the
%   parameters manually
%
% - If Use_scenario = 1 : we use an existing scenario and load the
%   parameters from a file.txt
%
% - If Use_scenario has a different value : we do not use an existing scenario, and 
%   we give default input parameters

Use_scenario = 0;

% We here choose the existing scenario we want to run (not necessary if Use_scenario ~= 1):
num_scenario = 13;



if Use_scenario == 1 
    %% WE RUN AN EXISTING SCENARIO:
    
    % We get the parameters from a .txt file, and save the results in anothe text file...

    % We open the file containing the parameters:
    dest = 'Experiences\Scenario_';
    dest = strcat(dest,num2str(num_scenario));
    dest = strcat(dest,'\');
    addpath(dest) %% So that we only add to the path the folder containing the scenario we want to run
    fileID = fopen( strcat(dest,'Parameters.txt') , 'r' );
    formatSpec = '%f';

    % We extract the parameters from the file (see below for the details on the parameters):
    param = fscanf(fileID,formatSpec);

    J          = param(1);
    n          = 2^J;
    J_c        = param(2);
    c          = 2^J_c;
    M          = param(3);
    V          = param(4);
    delta_p    = param(5);
    slopeRange = param(6);
    w          = param(7);
    ymin       = param(8);
    ymax       = n - param(9);
    csc        = param(10);
    fsc        = param(11);
    maxLength  = param(12); % Calculates BP for lengths 1,...,maxLength

    formatSpec = '%s';
    param = fscanf(fileID,formatSpec);
    fclose(fileID);
    type_image = param;

    
    % Parameters of the horizon curve :
    C_alpha = 2;
    alpha   = 1.7;
    


elseif Use_scenario == 0;    
  %% THEN WE DO NOT RUN AN EXISTING SCENARIO AND GIVE THE PARAMETERS MANUALLY:

    % We manually enter the parameters through the Main_EdgeletLab.m...   

    % Size of the image :
    J = 7;                                      % dyadic size of image
    n =  2^J;                                   % size of image (Note that the image is a square).
    
    % Sampling parameters :
    J_c = 7;                                    %(must be <= J)
    c   = 2^J_c;                                % number of columns we sample

    % Mean and variance of the additive gaussian noise:
    M = 0;                                      % mean
    V = 1;                                      % variance :

    % Vertical distance between two ends of edgelets on the same vertical line :
    delta_p = 1;                                % (in pixel) It must be of the form : delta_p = 1 / 2^k, with k is an integer
    delta = delta_p / n;                        % (in distance)

    % Maximum slope of an edgelet (in absolut value) :
    slopeRange = 1;                             % (45 degrees slope, if we consider the entire image)
    slopeRange = slopeRange * 2^(J - J_c);      % we adjust the slope considering the sampling 

    % depth of the filter :
    w = 5;                                      % (in pixels)

    % Upper bound and lower bound of the horizon function on the y-axis (in pixels) to not compute all the image :
    ymin = w;
    ymax = n - w;


    % Coarsest and finest scales :
    csc = J - J_c;                              % Coarsest scale we can compute (warning :must be >= J - J_c) 
%     csc = 2;
    fsc = J-1;                                  % finest scale we can compute 

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



    
else
    %% THEN WE DO NOT RUN AN EXISTING SCENARIO AND GIVE DEFAULT INPUT PARAMETERS:
    
    % The parameters have not been given to the algorithm, so we give default values... 
    J   = 5;
    n   = 2^J;
    M   = 0; 
    V   = 0.2;
    type_image = 'sin';
    maxLength = 5;
    disp('Error : The variable "Use_scenario" must be equal to 0 or 1 in Main_EdgeletLab.m');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE AN IMAGE
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

% Plot the noiseless image:
figure(1)
imshow(I_temp);

% Plot the noisy image:
% figure(2)
% imshow(I_b_temp);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Canny filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We perform a Canny filter :

thresh  = [0.2 , 0.6];                % Canny filter threshold
I_canny = edge(I_b,'canny',thresh);   % We perform the Canny filter

% [I_canny, thresh] = edge(I_b,'canny');   % We calculate the threshold
% I_canny = edge(I_b,'canny', thresh);   % We perform the Canny filter

figure;
imshow(I_canny);


% We highlight in red the pixels recognized by Canny :

figure;
rgbImage = highlight_red( I_b_temp, I_canny );
imshow(rgbImage);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Taking the parameters...');
if Use_scenario == 0 || Use_scenario == 1
    graphparam = GetEdgeletGraphParam(n, csc, fsc, delta_p, slopeRange, ymin, ymax, C_alpha, alpha, w, c);
else
    graphparam = GetEdgeletGraphParam(n);
end



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

% Threshold :
N_j = 2 * n^2 / delta_p^2 ./ [ 1 : maxLength ];
N_paths = N_j .* ( 2*n / delta_p ./ [ 1 : maxLength ]) .^ ([ 1 : maxLength ] - 1);
alpha_s = 1; % slowly depends on t_s ... So alpha_s = 1 is what we use in practice.
% t_s = 0;   % threshold for the best paths. We consider the best path only if its cost is > t_s
t_s = sqrt( ( V .* ([ 1 : maxLength ]) ) .* ( 2*log(N_paths) - log(log( N_paths )) - log(4*pi) - 2*log(alpha_s) ) );

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
%       figure();
%       imshow(binary_image);
      
      DisplayEdgeletPath(bestpaths{k},graphparam, sample); % plot the best path

      title(strcat('Best path restricted to #edgelets=',num2str(k),...
                       ', |cost|=',num2str(abs(costpath(k)))));
                   
        
      if Use_scenario == 1
      % Then we use an existing scenario...
          
      % Then we save the results :
        fileID = fopen(strcat(dest,'Results.txt'),'w');
        fprintf(fileID,'Best path : ');
        fprintf(fileID,'%f\n ',bestpaths{k});
        fclose(fileID);
        
      % And we save the figure k:   
        name = num2str(k);
        name = strcat('Path_',name,'_edgelets'); % we write the name of the figure k
        saveas(h, strcat(dest,name),'fig');
        
      end
      
  else
      fprintf( 'Warning : There is no path of length l = %d \n', k);
  end
  
end
 
disp('Done!');


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: Main_EdgeletLab.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego