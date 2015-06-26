% Gruppenmitglieder: Simon Bilgeri, Christoph Wittmann

%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'k',0.05,'N',100,'min_dist',80,'segment_length',9,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'k',0.05,'N',100,'min_dist',80,'segment_length',9,'do_plot',false);


%% KorrespondenzschÃ¤tzung
W=25;
min_corr=0.95;
tic
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,W,min_corr,'do_plot',false);
t=toc;

disp(['Es wurden ' num2str(size(Korrespondenzen,2)) ' Korrespondenzpunktpaare in ' num2str(t) 's gefunden.'])

%% Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus'
Korrespondenzen_robust = F_ransac(Korrespondenzen,'tolerance',0.015);


%% Zeige die robusten Korrespondenzpunktpaare
figure('name', 'Punkt-Korrespondenzen nach RANSAC');
imshow(uint8(IGray1))
hold on
plot(Korrespondenzen_robust(1,:),Korrespondenzen_robust(2,:),'r*')
imshow(uint8(IGray2))
alpha(0.5);
hold on
plot(Korrespondenzen_robust(3,:),Korrespondenzen_robust(4,:),'g*')
for i=1:size(Korrespondenzen_robust,2)
    hold on
    x_1 = [Korrespondenzen_robust(1,i), Korrespondenzen_robust(3,i)];
    x_2 = [Korrespondenzen_robust(2,i), Korrespondenzen_robust(4,i)];
    line(x_1,x_2);
end
hold off


%% Berechne die Essentielle Matrix
load('K.mat');
E = achtpunktalgorithmus(Korrespondenzen_robust,K);


%% Extraktion der möglichen euklidischen Bewegungen aus der Essentiellen Matrix und 3D-Rekonstruktion der Szene
 [T1,R1,T2,R2] = TR_aus_E(E);
 [T,R,lambdas,P1] = rekonstruktion(T1,T2,R1,R2,Korrespondenzen_robust,K);

%% Berechnung der Rückprojektion auf die jeweils andere Bildebene
% repro_error = rueckprojektion(Korrespondenzen_robust, P1, IGray2, T, R, K);


