function [Fx,Fy]=sobel_xy(Image)
% Performs sobel filter and return Fx, Fy
    Gray_image=rgb_to_gray(Image);
    Sobel_h=[1,0,-1;2,0,-2;1,0,-1];
    Sobel_v=[1,2,1;0,0,0;-1,-2,-1];
    % perform convolution, return same size as image 'same' and normalize
    Fx=(conv2(double(Gray_image),Sobel_h,'same')*(1/8)*log(2)); 
    Fy=(conv2(double(Gray_image),Sobel_v,'same')*(1/8)*log(2));
end

