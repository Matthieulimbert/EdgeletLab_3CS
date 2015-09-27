%
% Getting_threshold.m
% 
% This Matlab scripts gets the detection threshold to use in the Matlab file "Experiment_detection_rate.m".
% You must run EdgeletPath.m before running this file.
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
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cleaning the workplace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clear path

disp('--Experiment_detection_rate.m--');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION OF THE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

  %WE DO NOT RUN AN EXISTING SCENARIO AND GIVE THE PARAMETERS MANUALLY:

    % We manually enter the parameters through the Main_EdgeletLab.m...   

    % Size of the image :
    J = 9;                                      % dyadic size of image
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
%     csc = J - J_c;                                        % Coarsest scale we can compute (warning :must be >= J - J_c) 
% %     csc = 2;
%     fsc = J-1;                                            % finest scale we can compute
    number_edgelets_path = 8;                               % we are in the monoscale case
    csc = log(number_edgelets_path) / log(2) + J - J_c;     % Coarsest scale we can compute (warning :must be >= J - J_c) 
    fsc = csc;                                              % finest scale we can compute 



    % Type of image :
    type_image = 'noise_only';


    % Maximum size of the path :
    maxLength = 8;                              % Calculates BP for lengths 1,...,maxLength
    
    
    % Parameters of the horizon curve :
    C_alpha = 2;
    alpha   = 1.7;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STARTING OF THE LOOP ON THE VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pas = 0.5;
V_min = 0;
V_max = 2;
threshold = zeros( 1, 1 + (V_max - V_min) / pas );

indice_V = 1; % initialization of the indice on the variance


for V = V_min:pas:V_max
    
fprintf( 'Starting on the loop for V = %d \n', V);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STARTING OF THE LOOP ON THE ITERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_it = 3;
cost_vector = zeros ( 1, nb_it ); % initialization of the counter "cost"

for it = 1:nb_it

fprintf( '\nIteration %d / %d for V = %d \n', it, nb_it, V);
    
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WE UPDATE THE VECTOR CONTAINING THE COSTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cost_vector(it) = costpath(number_edgelets_path);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF THE LOOP ON THE ITERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATION OF THE THRESHOLD FOR A CERTAIN VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

rejection_pourcentage = 0.95; % pourcentage of rejection in case there is no path...
indice_threshold      = round(rejection_pourcentage * nb_it);
cost_vector = sort(cost_vector);
threshold(indice_V) = cost_vector(indice_threshold);
indice_V = indice_V + 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF THE LOOP ON THE VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We save the result of the experience
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Path of the file containting the results :
pourcentage = 2^J_c / 2^J * 100; % pourcentage of the image computed
dest = 'Experiences\Detection_rate\';
addpath(dest) %% So that we only add to the path the folder containing the scenario we want to run


% Name of the file containting the results :
name = 'Thresholds_';
name = strcat(name,num2str(pourcentage));
name = strcat(name,'.txt');

fileID = fopen(strcat(dest,name),'w');


% We write the thresholds that we computed :
for i = 1 : length(threshold)
fprintf(fileID,'%f\n', threshold(i));
end


% We write the parameters
fprintf(fileID,'\n\n\n PARAMETERS : \n\n');
fprintf(fileID,'image type : ');
fprintf(fileID,'%s\n', type_image);
fprintf(fileID,'image size J: ');
fprintf(fileID,'%d\n', J);
fprintf(fileID,'image size computed J_c: ');
fprintf(fileID,'%d\n', J_c);
fprintf(fileID,'number of edgelets in the path: ');
fprintf(fileID,'%d\n', number_edgelets_path);
fprintf(fileID,'number of iterations nb_it: ');
fprintf(fileID,'%d\n', nb_it);
fprintf(fileID,'filter window size w : ');
fprintf(fileID,'%d\n', w);
fprintf(fileID,'Variance vector : \n');
fprintf(fileID,'%f\n\n', V_min:pas:V_max);
fclose(fileID);


disp('Done!');


% Note: Code was adapted from ChirpLab Version 1.1

% $RCSfile: Main_EdgeletLab.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego