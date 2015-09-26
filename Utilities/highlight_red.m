function rgbImage = highlight_red( I_1, I_2 )

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


rgbImage = cat(3,I_1, I_1, I_1); % we convert the grayscale image into a rgb image


% We plot in red the pixels highlighted whose coordinates are in I_2 :

rgbImage1               = rgbImage(:,:,1);
rgbImage1(I_2 == 1) = 1;
rgbImage(:,:,1)         = rgbImage1;

rgbImage1               = rgbImage(:,:,2);
rgbImage1(I_2 == 1) = 0;
rgbImage(:,:,2)         = rgbImage1;

rgbImage1               = rgbImage(:,:,3);
rgbImage1(I_2 == 1) = 0;
rgbImage(:,:,3)         = rgbImage1;

