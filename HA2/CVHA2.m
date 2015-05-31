% Gruppenmitglieder: Simon Bilgeri, Christoph Wittmann

%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'k',0.05,'N',100,'min_dist',80,'segment_length',9,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'k',0.05,'N',100,'min_dist',80,'segment_length',9,'do_plot',false);


%% Korrespondenzschätzung
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2, 'window_length', 121, 'min_corr', 0.90, 'do_plot', true);