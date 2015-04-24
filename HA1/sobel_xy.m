function [Fx,Fy]=sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zur�ckgibt.
Gray_image=rgb_to_gray(Image);
Sobel_h=[1,0,-1;2,0,-2;1,0,-1];
Sobel_v=[1,2,1;0,0,0;-1,-2,-1];
Fx=(conv2(Gray_image,Sobel_h,'same')*(1/8)*log(2));
Fy=(conv2(Gray_image,Sobel_v,'same')*(1/8)*log(2));
end

