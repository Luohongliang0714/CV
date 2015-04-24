function [Gray_image]=rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es zurückgegeben werden.
if size(Image,3)==1
    Gray_image=Image;
else
    Image=double(Image); % Konvertierung für double Genauigkeit
    % Jeder Kanal wird mit Faktor multipliziert und aufaddiert
    Gray_image=uint8(Image(:,:,1)*0.299+Image(:,:,2)*0.587+Image(:,:,3)*0.114);
end
