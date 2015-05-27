function [Gray_image]=rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es zurückgegeben werden.

% Überprüfe, ob das Bild tatsächlich drei Kanäle hat
if(size(Image,3)~=3)
    Gray_image=Image;
    return;
end
% Umwandlung nach ITU-Formel
Gray_image = 0.299*Image(:,:,1)+0.587*Image(:,:,2)+0.114*Image(:,:,3);

end
