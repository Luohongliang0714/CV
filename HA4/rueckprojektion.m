function repro_error=rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler

%% Reproject 3D points onto second camera

% Transformiere Korrespondenzen in Bild 2 in homogene Koordinaten:
x2 = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
%x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
% Construct R|T Matrix:
R_T=[R,T];
% P1 in homogeneous coordinates:
P1_hom=[P1;ones(1,size(Korrespondenzen,2))];
% Backprojection
x2_back=K*R_T*P1_hom;
% normalization
x2_back=x2_back./(ones(3,1)*x2_back(3,:));
% mean reprojection error
repro_error=sum(sqrt(sum((x2_back-x2).^2,1)))/size(x2,2);
%x1_back=K*P1_hom(1:3,:);
%x1_back=x1_back./(ones(3,1)*x1_back(3,:));

display(repro_error);

show_back_proj(I2,x2,x2_back);

end