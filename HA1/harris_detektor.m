function  Merkmale=harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

 %% Input parser
    P = inputParser;

    % Liste der notwendigen Parameter
    % Ein Bild als Input ist zwingend notwendig
    P.addRequired('Image', @(x)isnumeric(x));

    % Liste der optionalen Parameter
    % Größe des Bildsegments
    % (default) = 10;
    P.addOptional('segment_length', 10, @(x) x ~=0);
    % Ecken- und Kantengewichtung mit k
    % (default) = 0.05;
    P.addOptional('k', 0.05, @(x)isscalar(x)&isnumeric(x));
    % Schwellenwert tau
    % (default) = 1;
    P.addOptional('tau', 0.5, @(x)isscalar(x)&isnumeric(x));
    P.addOptional('do_plot', 1, @(x)islogical(x));
    P.addOptional('tile_size', 10, @(x) (isscalar(x) | prod(size(x)==[2,1]))&isnumeric(x));
    P.addOptional('N', 20, @(x)isscalar(x)&isnumeric(x));
    P.addOptional('min_dist', 30, @(x)isscalar(x)&isnumeric(x));
    % Lese den Input
    P.parse(Image,varargin{:});

    %% Extrahiere die Variablen aus dem Input-Parser
    Img     = P.Results.Image; 
    % Img     = Bild;    % Auch moeglich, da die Variable Bild der Funtion direkt uebergeben wurde  
    slength   = P.Results.segment_length;
    k   = P.Results.k;
    tau = P.Results.tau;
    do_plot = P.Results.do_plot;
    tile_size=P.Results.tile_size;
    N=P.Results.N;
    min_dist=P.Results.min_dist;
    sigma=1;
    %% Harris detector
    %% Convert to greyscale and perform sobel filter
    Img_gray=rgb_to_gray(Img);
    [Fx,Fy]=sobel_xy(Img_gray);
    %% Create Covariance matrix for each pixel to approximate Harris matrix
    G=zeros(size(Img,1),size(Img,2),4);
    H=zeros(size(Img,1),size(Img,2));
    H_mask=zeros(size(Img,1),size(Img,2));
    G(:,:,1)=Fx.^2;
    G(:,:,2)=Fx.*Fy;
    G(:,:,3)=G(:,:,2);
    G(:,:,4)=Fy.^2;
    segment_half = floor(slength/2);            % Halbe Fesntergroesse
    ind = -segment_half:segment_half;               % indices (Bsp.: -1,0,1)
    gau = exp(-(ind.^2)/(2*sigma^2))';    % 1D Gauss
    Weight_filter = gau;                           % 2 Separable Filter ergeben das 2D Filter
    Weight_filter_1d = Weight_filter./sum(Weight_filter);                 % Summe aller Koeffizienten ergibt 1
    for i=1:4
        G(:,:,i) = conv2(double(G(:,:,i)),Weight_filter_1d,'same');
        G(:,:,i) = conv2(double(G(:,:,i)),Weight_filter_1d','same');
    end
    % Calculate H Value for each pixel matrix, first determinant, then
    % trace
    H=G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3)-k*(G(:,:,1)+G(:,:,4)).^2;
    %H=(H>-tau) .* (H < tau); %Textureless
    %H=H<-tau; %Kante
    
    H_mask=H>tau; % Ecke

    % free memory
    clearvars G;
    % Filter H matrix with mask
    H=H.*H_mask;
    % reshape H matrix with tile_size
%     if isscalar(tile_size)
%         randy=mod(size(H,1),tile_size);
%         randx=mod(size(H,2),tile_size);
%         num_tiles=((size(H,1)-randy)/tile_size)*((size(H,2)-randx)/tile_size);
%         Tiles=reshape(H(1:size(H,1)-randy,1:size(H,2)-randx),tile_size,tile_size,num_tiles);
%     else
%         % falls vector
%         Tiles=reshape(H,tile_size(1),tile_size(2));
%     end
    % Filter all Features which do not have the min_dist. 
    % Search for next_neighbor and keep larget one
    [feat_r,feat_c]=find(H);
    feat_mat=[feat_r,feat_c]';
    feat_mat_res=feat_mat;
    dist=zeros(size(feat_mat));
    dist2=zeros(size(feat_mat,2));
    for i=1:size(feat_r,1)
        dist(1,:)=(feat_mat(1,:)-feat_mat(1,i))^2;
        dist(2,:)=(feat_mat(2,:)-feat_mat(2,i))^2;
        dist2=dist(1,:)+dist(2,:); % Wurzel brauchen wir nicht
        dist2(dist2==0)=[];
        nn_ind=min(dist2)
    end
    num_tilesy=((size(H,1)-randy)/tile_size);
    num_tilesx=((size(H,2)-randx)/tile_size);
    M=zeros(2,N*num_tilesy*num_tilesx);
    corner_count=0;
    for i=1:num_tilesy
        for j=1:num_tilesx
            tile=H(1+(i-1)*tile_size:1+i*tile_size,1+(j-1)*tile_size:1+j*tile_size);
            [sortedValues,sortIndex] = sort(tile(:),'descend');
            maxIndex = sortIndex(1:N);
            maxIndex(sortedValues(1:N)==0)=[];
            [row,col] = ind2sub(size(tile), maxIndex);
            row=row+(i-1)*tile_size;
            col=col+(j-1)*tile_size;
            M(:,corner_count+1:corner_count+size(row,1))=[col,row]';
            corner_count=corner_count+size(row,1);
        end
    end
    M=M(:,any(M));
    Merkmale=M;

        %% show corners if specified
    if do_plot
        figure;
        title('Harris Detector using tiles and min distance');
        imshow(Img)
        hold on;            %# Add subsequent plots to the image
        plot(M(1,:),M(2,:),'rx','MarkerSize',5);  %# NOTE: x_p and y_p are switched (see note below)!
        hold off; 
        
    end
end

