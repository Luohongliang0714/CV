% Gruppenmitglieder: Christoph Wittmann, Simon Bilgeri

%F�r die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%enfernen und sicherstellen, dass alle optionalen Parameter über den
%entsprechenden Funktionsaufruf fun('var',value) modifiziert werden können.

%% Bild laden
 Image = imread('szene.jpg');
 IGray = rgb_to_gray(Image);


%% Harris-Merkmale berechnen
 tic
 %Merkmale = harris_detektor(IGray,'do_plot',true);
 Merkmale = harris_detektor(Image,'tau',5000,'do_plot',true,'tile_size',36,'N',3,'segment_length',3);
 toc
