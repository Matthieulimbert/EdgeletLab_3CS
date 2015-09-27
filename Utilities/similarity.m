function dist = similarity( I_1, I_2 )

% highlight_red --  This function highlights some pixels (by plotting them in red)
% of the matrix I_1 whose coordinates are in I_2.
% I_1 is a grayscale image, and I_2 is a binary image with I_2(x,y) = 1 if we
% want to highlight the pixel (x,y), and I_2(x,y) = 0 if we don't
% want to highlight the pixel (x,y).
%
%  Usage
%   rgb = highlight_red( I_1, I_2 )
%  Inputs
%    I_1		 graycale image  
%    I_2         binary matrix containing the coordinates of the pixels of
%    I_1         to highlight. I_2(x,y) = 1 if we want to highlight the pixel 
%                (x,y), and I_2(x,y) = 0 if we don't want to highlight the pixel (x,y).
%  Outputs
%    rgbImage    rgb image. This is the grayscale image I_1 with the pixels
%                highlited in red.
%


n = length(I_1); % size of the image

coord = (1:n)';
coord = repmat(coord, 1, n); % matrix containing the vertical coordinates

I_1   = I_1  (:,2:n-1);
I_2   = I_2  (:,2:n-1);
coord = coord(:,2:n-1);

I_1 = I_1 .* coord;
I_2 = I_2 .* coord;

I_1 = max(I_1);
I_2 = max(I_2);

dist = max( abs(I_1 - I_2) );



