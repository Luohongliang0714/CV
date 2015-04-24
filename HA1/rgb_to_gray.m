function [Gray_image]=rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es zurückgegeben werden.
if size(Image,3)==1
    Gray_image=Image;
else
    % Way too slow!
%     Image=double(Image);
%     Gray_image=zeros(size(Image,1),size(Image,2));
%     for x=1:size(Image,1)
%         for y=1:size(Image,2)
%             Gray_image=squeeze(Image(x,y,:))*[0.299,0.587,0.114];
%         end
%         fprintf('%d',x);
%     end
    % Multiply each channel with its factor and sum them together
    Image=double(Image);
    Gray_image=uint8(Image(:,:,1)*0.299+Image(:,:,2)*0.587+Image(:,:,3)*0.114);

end
