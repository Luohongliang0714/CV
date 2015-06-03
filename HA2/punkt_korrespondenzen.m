function [ Korrespondenzen ] = punkt_korrespondenzen( I1, I2, Mpt1, Mpt2, varargin) 
%function: punkt_korrespondenzen
%INPUT: 2 Bilder I1,I2; Feature Punkte Mpt1, Mpt2, und optionale Parameter.
%OUTPUT: Korrespondenzmatrix in Form [x_p1,y_p1,x_p2,y_p2]

%% Input parser
P = inputParser;

% Liste der optionalen Parameter
% Fensterlänge
P.addOptional('window_length', 5, @isnumeric)
% Minimum Correlations-Wert für Match
P.addOptional('min_corr', 0.8, @isnumeric)
% Plot ein/aus
P.addOptional('do_plot', false, @islogical);

% Lese den Input
P.parse(varargin{:});

% Extrahiere die Variablen aus dem Input-Parser
window_length = P.Results.window_length;
min_corr = P.Results.min_corr;
do_plot = P.Results.do_plot;
% Check if window_length is odd, if not change it to an odd number
if mod(window_length,2)==0
    window_length=window_length+1;
    fprintf('Der Parameter window_length sollte eine ungerade Ganzzahl sein, nutze deshalb window_length = %i\n',window_length);
end
%% NCC Matching
tic
num_Merkmale1=size(Mpt1,2);
num_Merkmale2=size(Mpt2,2);
Korrespondenzen=[];
window_size=floor(window_length/2);
% Extract intesity values in window for all features in both images
intensity_val1=-1*ones(window_length^2,num_Merkmale1);
intensity_val2=-1*ones(window_length^2,num_Merkmale2);
num_merkmal_border=[0;0];
for i=1:num_Merkmale1
    try
        intensity_val1(:,i)=reshape(I1(Mpt1(2,i)-window_size:Mpt1(2,i)+window_size,Mpt1(1,i)-window_size:Mpt1(1,i)+window_size),window_length^2,1);
    catch
        num_merkmal_border(1)=num_merkmal_border(1)+1;
    end
end
for j=1:num_Merkmale2
    try 
        intensity_val2(:,j)=reshape(I2(Mpt2(2,j)-window_size:Mpt2(2,j)+window_size,Mpt2(1,j)-window_size:Mpt2(1,j)+window_size),window_length^2,1);
    catch
        num_merkmal_border(2)=num_merkmal_border(2)+1;
    end
end
% Discard Features too close to the borders (-1 intensity values)
Mpt1(:,(intensity_val1(1,:)==-1))=[];
intensity_val1(:,(intensity_val1(1,:)==-1))=[];
Mpt2(:,(intensity_val2(1,:)==-1))=[];
intensity_val2(:,(intensity_val2(1,:)==-1))=[];
% Precompute mean and variances:
mean_val1=mean(intensity_val1);
mean_val2=mean(intensity_val2);
stdev_val1=sqrt(var(intensity_val1));
stdev_val2=sqrt(var(intensity_val2));
num_Merkmale1=size(Mpt1,2);
% normalize intensity values
intensity_val1=bsxfun(@rdivide,bsxfun(@minus,intensity_val1,bsxfun(@times,ones(window_length^2,1),mean_val1)),bsxfun(@times,ones(window_length^2,1),stdev_val1));
intensity_val2=bsxfun(@rdivide,bsxfun(@minus,intensity_val2,bsxfun(@times,ones(window_length^2,1),mean_val2)),bsxfun(@times,ones(window_length^2,1),stdev_val2));
NCC=1/((window_length)^2-1)*intensity_val1'*intensity_val2;
[ncc_val,ncc_id]=max(NCC,[],2);
ncc_id_m2=ncc_id(ncc_val>=min_corr);
ncc_id_m1=(1:num_Merkmale1)';
ncc_id_m1=ncc_id_m1(ncc_val>=min_corr);
Korrespondenzen=[Mpt1(:,ncc_id_m1);Mpt2(:,ncc_id_m2)];
t1=toc;

if (sum(num_merkmal_border)>0)
    fprintf('Achtung: %i Merkmale in Bild I und %i in Bild II liegen zu dicht am Rand und werden nicht berücksichtigt.\n', num_merkmal_border);
end
fprintf('%i Korrespondenzen wurden in %.3g Sekunden gefunden\n', size(Korrespondenzen,2),t1);

%% Visualize found matches
if do_plot
    % Sortiere Korrespondenzen nach x1 und y1:
    Korrespondenzen=sortrows(Korrespondenzen',[1,2])';
    colors=char('blue','red','green','magenta','yellow','cyan');
    markers=char('o','s','^','v','>','<','p','h');
    marker_size=7;
    font_size=8;
    offset=2*font_size+marker_size+1;
    figure;
    imshow([I1,I2]);
    hold on;
    title('Found matches')
    for i=1:size(Korrespondenzen,2)
        plot(Korrespondenzen(1,i),Korrespondenzen(2,i),strcat(colors(mod(i,size(colors,1))+1,:),markers(mod(i,size(markers,1))+1,:)),'MarkerSize',marker_size)
        text(Korrespondenzen(1,i),Korrespondenzen(2,i)+offset,num2str(i),'Color',colors(mod(i,size(colors,1))+1,:),'FontSize',font_size)
        plot(Korrespondenzen(3,i)+size(I1,2),Korrespondenzen(4,i),strcat(colors(mod(i,size(colors,1))+1,:),markers(mod(i,size(markers,1))+1,:)),'MarkerSize',marker_size)
        text(Korrespondenzen(3,i)+size(I1,2),Korrespondenzen(4,i)+offset,num2str(i),'Color',colors(mod(i,size(colors,1))+1,:),'FontSize',font_size)
        line([Korrespondenzen(1,i),Korrespondenzen(3,i)+size(I1,2)],Korrespondenzen([2,4],i),'Color',colors(mod(i,size(colors,1))+1,:));
    end
end
end

