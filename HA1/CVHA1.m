% Gruppenmitglieder: Christoph Wittmann, Simon Bilgeri

% Infos:
% Die Detektor Werte sind mit dem größten Wert normiert, Tau also relativ.
% (Tau_max=1000)
% Extra Funktion find_N_features.m ber
%% Bild laden
Image = imread('szene.jpg');
IGray = rgb_to_gray(Image);

%% Harris-Merkmale berechnen
tic
Merkmale = harris_detektor(IGray,'k',0.05,'tau',1,'do_plot',true,'tile_size',[400;400],'N',3,'segment_length',4,'min_dist',5);
toc
fprintf('#Anzahl der Merkmale: %d\n', size(Merkmale,2));
 
%% Checkerboard test
C=checkerboard(30);
tic
Merkmale_C=harris_detektor(C,'k',0.05,'tau',1,'do_plot',true,'tile_size',10,'N',3,'segment_length',4,'min_dist',5);
toc
fprintf('#Anzahl der Schachbrettecken: %d\n', size(Merkmale_C,2));