function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1

% Transformiere Korrespondenzen in homogene und kalibrierte Koordinaten:
x1 = K\[Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
x2 = K\[Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];

T_mat={T1,T2};
R_mat={R1,R2};

% Stelle für jede der 4 Möglichkeiten Matrix M auf für 1. und 2. Bild und
% zähle die Anzahl der positiven lambdas:
max_pos_lambdas=0;
for r=1:2
    for t=1:2
        % M11
        % Bild 1
        M=cross(x2,R_mat{r}*x1);
        cell_vec=mat2cell(M,3,ones(size(Korrespondenzen,2),1));
        M=[blkdiag(cell_vec{:}),reshape(cross(x2,T_mat{t}*ones(1,size(Korrespondenzen,2))),size(Korrespondenzen,2)*3,1)];
        [~,~,V]=svd(M);
        lambdas_temp=V(1:end-1,end);
        % Bild 2
        M=cross(x1,R_mat{r}'*x2);
        cell_vec=mat2cell(M,3,ones(size(Korrespondenzen,2),1));
        M=[blkdiag(cell_vec{:}),reshape(cross(x1,R_mat{r}'*(-1)*T_mat{t}*ones(1,size(Korrespondenzen,2))),size(Korrespondenzen,2)*3,1)];
        [~,~,V]=svd(M);
        lambdas_temp=[lambdas_temp,V(1:end-1,end)];
        % num_positive lambdas:
        if (sum(lambdas_temp(:)>0) > max_pos_lambdas)
            max_pos_lambdas=sum(lambdas_temp(:)>0);
            lambdas=lambdas_temp;
            R=R_mat{r};
            T=T_mat{t};
        end
    end
end

P1=(ones(3,1)*lambdas(:,1)').*x1;

% Plot resulting camera frames and points
plot3(P1(1,:),P1(2,:),P1(2,:),'o');
end