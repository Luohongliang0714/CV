% Gruppenmitglieder: Christoph Wittmann, Simon Bilgeri

%Für die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%enfernen und sicherstellen, dass alle optionalen Parameter Über den
%entsprechenden Funktionsaufruf fun('var',value) modifiziert werden können.

%% Bild laden
 Image = imread('szene.jpg');
 IGray = rgb_to_gray(Image);


%% Harris-Merkmale berechnen
 tic
 %Merkmale = harris_detektor(IGray,'do_plot',true);
 Merkmale = harris_detektor(IGray,'tau',1000,'do_plot',true,'tile_size',[400;400],'N',3,'segment_length',4,'min_dist',5);
 toc
 fprintf('#Anzahl der Merkmale: %d\n', size(Merkmale,2));
 
%% Checkerboard test
C=checkerboard(20);
tic
Merkmale_C=harris_detektor(C,'tau',1000,'do_plot',true,'tile_size',10,'N',3,'segment_length',4,'min_dist',5);
toc
fprintf('#Anzahl der Schachbrettecken: %d\n', size(Merkmale_C,2));