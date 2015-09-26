function [ im ] = load_image( name_image, n )
% DisplayChirplets -- Converts the image in a matrix corresponding of the
%                     black and white image of size n x n
%            
%  Usage
%    [ im ] = load_image( name_image, n )
%  Inputs
%    name_image       name of the image we want to load
%    n                size of the image we want to compute (NB the image
%                     can be bigger than the size n x n. In such case, the
%                     image is cropped). The image loaded must be of size n
%                     x n or bigger
%  Outputs
%    im               The matrix corresponding to the n x n black and white image
% 
%    
%  See Also
%    - MakeImage



% We read the image:
im = imread(name_image);
im = rgb2gray(im); % RGB image converted in black and white image
im_1 = mat2gray(im); % Values of the matrix are now between 0 (black) and 1 (white).
im = 1 - im_1; % Now pixel = 1 means it is black, and pixel = 0 means it is white.


% We crop the image:
im   = im (1 : n, 1 : n);

% We display the image :
% im_1 = im_1 (1 : n, 1 : n);
% imshow(im_1);

end


% $RCSfile: load_image.m,v $
% $Date: 06/23/2015 $
% $Revision: 1 $
%
% Copyright (c) Matthieu Limbert, University of California, San Diego