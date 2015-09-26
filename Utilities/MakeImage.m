function [I, I_b] = MakeImage(n, M, V, type_image)

% MakeChirp -- Make artificial image
%
%   Usage
%     [I, I_b] = MakeImage(n, M, V, type_image)
%   Inputs
%     n             desired image size
%     M             mean of the gaussian additive noise
%     V             variance of the gaussian additive noise
%     type_image    it is either :
%                   - a name of an image that exists in the path (and then
%                   we load it). example : 'Road2.bmp'.
%                   - a type of image ('line', 'sin' or 'sin_am') that we
%                   directly build in this function
%   Outputs
%     I             noiseless black and white image (matrix n x n)
%     I_b           noisy black and white image (matrix n x n)
%
%  See Also
%    - load_image


B = 1;
W = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We build the noiseless image, depending on the type:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = zeros(n,n);




if strcmp(type_image, 'noise_only')
%% the image is constant, with only noise...


I = B * ones(n,n); % image is black

elseif strcmp(type_image, 'line')
%% The horizon function is a horizontal line on the middle of the image...

I(1:n/2,:) = 1;

I(I==1) = B;
I(I==0) = W;



elseif strcmp(type_image, 'sin') 
%% The horizon function is a horizontal curve on the middle of the image...
    
% Parameters to manage the sinus (horizon function) :

delta      = 0.4;
coeff      = 0.9;
coeff_mult = 0.7;

% We build the image :
f = 1 : 1 : n;
f_sin = sin ( f / (coeff * n) + delta ) * coeff_mult ;
f_sin = repmat(f_sin, n, 1);

coord = (1 : 1 : n)';
coord = (1/n) * repmat(coord, 1, n);

bool = coord - f_sin < 0;

I(bool == 1) = B;
I(bool == 0) = W;


elseif strcmp(type_image, 'sin_am')
%% The horizon function is a horizontal curve on the middle of the image...
    
% Parameters to manage the sinus (horizon function) :
T1 = n/2;
f0 = 0.05 / 3;

% We build the image :
t = 1 : 1 : n;
% f_sin = exp(-t/T1).*sin(2*pi*f0 * t);
f_sin = exp(-t/T1).* sin(2*pi*f0 * t);
f_sin = f_sin * 0.15 + 0.5;

f_sin = repmat(f_sin, n, 1);

coord = (1 : 1 : n)';
coord = (1/n) * repmat(coord, 1, n);

bool = coord - f_sin < 0;


I(bool == 1) = B;
I(bool == 0) = W;




else 
    %% We load an image that already exists:
    
    I = load_image( type_image, n );

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We build the corresponding noisy image:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise = normrnd(M,V,n,n); % gaussian distribution
I_b = I + noise; % additive gaussian noise

% I_b = imnoise(I,'gaussian',M,V);


end


% $RCSfile: load_image.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego