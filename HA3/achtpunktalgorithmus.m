function [EF]=achtpunktalgorithmus(Korrespondenzen,varargin)

 %% Input parser
    P = inputParser;
    
    P.addOptional('K', 0, @(x) size(x,1)==3 && size(x,2)==3);
    % (default) = 0; if 0 nicht gegeben
    %Kamerakalibrierungsmatrix

    P.parse(varargin{:});

    %% Extrahiere die Variablen aus dem Input-Parser
    K = P.Results.K; 
    if K==0
        compute_E=false;
    else
        compute_E=true;
    end


%% Korrespondenzen in homogenen Koordinaten
korr_hom1=Korrespondenzen(1:2,:);
korr_hom1(3,:)=ones(length(Korrespondenzen),1);

korr_hom2=Korrespondenzen(3:4,:);
korr_hom2(3,:)=ones(length(Korrespondenzen),1);

% Wenn K gegeben, kalibriere Koordinaten
if (compute_E)
    korr_hom1=K\korr_hom1;
    korr_hom2=K\korr_hom2;
end

%A-Matrix ( 4 Multiplikationen pro Korrespondenz )

A=[[korr_hom1(1,:);korr_hom1(1,:)].*korr_hom2(1:2,:);korr_hom1(1,:);[korr_hom1(2,:);korr_hom1(2,:)].*korr_hom2(1:2,:);korr_hom1(2,:);korr_hom2]';

%Singulärwertzerlegung von A Matrix
[U,S,V]=svd(A);

% Erzeuge Vektor aus 9ter Spalte, des V Vektors des kleinsten Singulärwerts
G=reshape(V(:,9),3,3);

% Erzeuge G Matrix
[U_g,S_g,V_g]=svd(G);

% Projeziere auf die Menge der E oder F Matrizen
if (compute_E)
    EF=U_g*[1,0,0;
            0,1,0;
            0,0,0]*V_g';
else
    EF=U_g*[S_g(1,1),0,0;
            0,S_g(2,2),0;
            0,0,0]*V_g';
end
  
end