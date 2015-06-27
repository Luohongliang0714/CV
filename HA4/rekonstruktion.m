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
        lambdas_temp=V(:,end);
        % Bild 2
        M=cross(x1,R_mat{r}'*x2);
        cell_vec=mat2cell(M,3,ones(size(Korrespondenzen,2),1));
        M=[blkdiag(cell_vec{:}),reshape(cross(x1,R_mat{r}'*T_mat{t}*ones(1,size(Korrespondenzen,2))),size(Korrespondenzen,2)*3,1)];
        [~,~,V]=svd(M);
        lambdas_temp=[lambdas_temp,V(:,end)];
        % num_positive lambdas:
        if (sum(sum(lambdas_temp(1:end-1,:)>0)) > max_pos_lambdas)
            max_pos_lambdas=sum(sum(lambdas_temp(1:end-1,:)>0));
            lambdas=lambdas_temp(1:end-1,:);
            R=R_mat{r};
            T=lambdas_temp(end,1)*T_mat{t};
        end
    end
end

% reconstruction
P1=(ones(3,1)*lambdas(:,1)').*x1;
P2=(ones(3,1)*lambdas(:,2)').*x2;
% Plot resulting camera frames and points
figure;

camera_length=0.005; %(0.5cm)
p1=[-camera_length;-camera_length;0];
p2=[-camera_length;+camera_length;0];
p3=[+camera_length;+camera_length;0];
p4=[+camera_length;-camera_length;0];
plot3([p1(3),p2(3),p3(3),p4(3),p1(3)],[p1(1),p2(1),p3(1),p4(1),p1(1)],[p1(2),p2(2),p3(2),p4(2),p1(2)],'b');
set(gca,'ydir','reverse')
hold on;
p1=R*p1 + T;
p2=R*p2 + T;
p3=R*p3 + T;
p4=R*p4 + T;
plot3([p1(3),p2(3),p3(3),p4(3),p1(3)],[p1(1),p2(1),p3(1),p4(1),p1(1)],[p1(2),p2(2),p3(2),p4(2),p1(2)],'r');
xlabel('Z');
ylabel('X');
zlabel('Y');
% draw cameras
plot3(0,0,0,'bo');
plot3(T(3),T(1),T(2),'ro');
plot3(P1(3,:),P1(1,:),P1(2,:),'gx');
end