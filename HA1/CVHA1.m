% Gruppenmitglieder: Christoph Wittmann, Simon Bilgeri

%F�r die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%enfernen und sicherstellen, dass alle optionalen Parameter �ber den
%entsprechenden Funktionsaufruf fun('var',value) modifiziert werden k�nnen.

%% Bild laden
 Image = imread('szene.jpg');
 IGray = rgb_to_gray(Image);


%% Harris-Merkmale berechnen
 tic
 %Merkmale = harris_detektor(IGray,'do_plot',true);
 Merkmale = harris_detektor(Image,'tau',1000,'do_plot',true,'tile_size',[400;400],'N',3,'segment_length',4,'min_dist',5);
 toc
 fprintf('#Anzahl der Merkmale: %d\n', size(Merkmale,2));
