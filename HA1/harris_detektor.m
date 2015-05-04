function  Merkmale=harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

 %% Input parser
    P = inputParser;

    P.addRequired('Image', @(x)isnumeric(x)); %Bild
    % Liste der optionalen ParameterMerkmale
    % Größe des Bildsegments
    % (default) = 4;
    P.addOptional('segment_length', 4, @(x) x ~=0);
    % Ecken- und Kantengewichtung mit k
    % (default) = 0.05;
    P.addOptional('k', 0.05, @(x)isscalar(x)&isnumeric(x));
    % Schwellenwert tau
    % (default) = 1;
    P.addOptional('tau', 1, @(x)isscalar(x)&isnumeric(x));
    P.addOptional('do_plot', 1, @(x)islogical(x));
    % Größe der Kacheln, default[ 200,200]
    P.addOptional('tile_size', 200, @(x) (isscalar(x) | prod(size(x)==[2,1]))&isnumeric(x));
    % Anzahl der Features pro Kachel N, default 3
    P.addOptional('N', 3, @(x)isscalar(x)&isnumeric(x)&x>1);
    % Minimaler Abstand der Features, default 3
    P.addOptional('min_dist', 3, @(x)isscalar(x)&isnumeric(x));
    % Lese den Input
    P.parse(Image,varargin{:});

    %% Extrahiere die Variablen aus dem Input-Parser
    Img = P.Results.Image; 
    slength = P.Results.segment_length;
    k = P.Results.k;
    tau = P.Results.tau;
    do_plot = P.Results.do_plot;
    tile_size=P.Results.tile_size;
    if isscalar(tile_size)
        tile_size=[tile_size;tile_size]; % mache skalar zu vector
    end
    N=P.Results.N;
    min_dist=P.Results.min_dist;
    sigma=1; % set to fixed value
    
    
    %% Harris detector
    %% Convert to greyscale and perform sobel filter
    Img_gray=rgb_to_gray(Img);
    [Fx,Fy]=sobel_xy(Img_gray);
    %% Create Jacobian matrix for each pixel to approximate Harris matrix
    fprintf('Calculate Harris Detector...');
    G=zeros(size(Img,1),size(Img,2),4);
    G(:,:,1)=Fx.^2;
    G(:,:,2)=Fx.*Fy;
    G(:,:,3)=G(:,:,2);
    G(:,:,4)=Fy.^2;
    segment_half = floor(slength/2);            % Halbe Fesntergroesse
    ind = -segment_half:segment_half;            % indices (Bsp.: -1,0,1)
    gau = exp(-(ind.^2)/(2*sigma^2))';    % 1D Gauss
    Weight_filter = gau;    % 2 Separable Filter ergeben das 2D Filter
    Weight_filter_1d = Weight_filter./sum(Weight_filter);      % Summe aller Koeffizienten ergibt 1
    G(:,:,1) = conv2(double(G(:,:,1)),Weight_filter_1d,'same');
    G(:,:,1) = conv2(double(G(:,:,1)),Weight_filter_1d','same');
    G(:,:,2) = conv2(double(G(:,:,2)),Weight_filter_1d,'same');
    G(:,:,2) = conv2(double(G(:,:,2)),Weight_filter_1d','same');
    G(:,:,3) = G(:,:,2); % symmetrische Matrix
    G(:,:,4) = conv2(double(G(:,:,4)),Weight_filter_1d,'same');
    G(:,:,4) = conv2(double(G(:,:,4)),Weight_filter_1d','same');
    
    % Calculate H Value for each pixel, first determinant, then
    % trace
    H=G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3)-k*(G(:,:,1)+G(:,:,4)).^2;
    % normalize H with largest value (comparable Tau)
    H=(1000/max(max(H)))*H;
    %H=(H>-tau) .* (H < tau); %Textureless
    %H=H<-tau; %Kante
    H_mask=H>tau; % Ecke
    fprintf('done.\n');
    % free memory, G not needed anymore
    clearvars G;
    
    %% Purge all smaller features within min_dist
    % Delate image with circular [min_dist*2+1,min_dist*2+1] mask
    % --> search for max element within window
    % only allow unchanged values -> smaller ones are purged
    % purge features at the border within min_dist first
    fprintf('Purge all smaller features in min_dist...');
    H_border=padarray(ones(size(H,1)-min_dist*2,size(H,2)-min_dist*2),[min_dist,min_dist]);
    H=H.*H_border;
    % Shift image (2*min_dist+1)^2-1 times
    shift_mat=zeros(size(H,1),size(H,2),2);
    shift_mat(:,:,1)=H;
    for i=-min_dist:min_dist
        for j=-min_dist:min_dist
            distance=sqrt(i^2+j^2);
            if ((j~=0 || i~=0) && distance <=min_dist)
                % shift matrix (borders dont care because they are 0 anyways)
                shift_mat(:,:,2)=circshift(H,[i,j]); % faster
                shift_mat(:,:,1)=max(shift_mat,[],3);
                % if the H values are the same, delete one of the features,
                H=H.*(shift_mat(:,:,2)~=H);
            end
        end
    end
    % Maximum where value did not change
    H_dil=(H==shift_mat(:,:,1));
    % Apply the filters
    H=H.*H_mask.*H_dil;
    fprintf('done.\n');
    %% Segment the image in tiles and take the N best features
    % Check if the tilesize is well chosen, if not only full tiles are used
    % first
    fprintf('Find N best features in each tile...');
    randy=mod(size(H,1),tile_size(2)); 
    randx=mod(size(H,2),tile_size(1));
    % number of uill tiles in x and y direction
    num_tilesy=((size(H,1)-randy)/tile_size(2));
    num_tilesx=((size(H,2)-randx)/tile_size(1));
    M=zeros(2,N*num_tilesy*num_tilesx); % preallocate for speed
    corner_count=0;
    for i=1:num_tilesy
        for j=1:num_tilesx
            % get the tile
            tile=H(1+(i-1)*tile_size(2):1+i*tile_size(2)-1,1+(j-1)*tile_size(1):1+j*tile_size(1)-1);
            [M,corner_count]=find_N_features(tile,N,M,i,j,corner_count,tile_size);
        end
    end
    

    %% Randbehandlung
    if randx~=0 
        %determine Nrx according to tile_size
        Nrx=ceil(N*((randx*tile_size(2))/(tile_size(1)*tile_size(2))));
        for i=1:num_tilesy
            %get the tile
            tile=H(1+(i-1)*tile_size(2):1+i*tile_size(2)-1,1+tile_size(1)*num_tilesx:tile_size(1)*num_tilesx+randx);
            [M,corner_count]=find_N_features(tile,Nrx,M,i,num_tilesx+1,corner_count,tile_size);
        end
    end
    
    if randy~=0
        %determine N according to tile_size
        Nry=ceil(N*((randy*tile_size(1))/(tile_size(1)*tile_size(2))));
        for j=1:num_tilesx
            %get the tile
            tile=H(1+tile_size(2)*num_tilesy:tile_size(2)*num_tilesy+randy, 1+(j-1)*tile_size(1):1+j*tile_size(1)-1);
            [M,corner_count]=find_N_features(tile,Nry,M,num_tilesy+1,j,corner_count,tile_size);
        end
    end
    
    % Check the last tile in the right corner
    Nrxy=ceil(N*((randx*randy)/(tile_size(1)*tile_size(2))));
    tile=H(1+tile_size(2)*num_tilesy:tile_size(2)*num_tilesy+randy,1+tile_size(1)*num_tilesx:tile_size(1)*num_tilesx+randx);
    M=find_N_features(tile,Nrxy,M,num_tilesy+1,num_tilesx+1,corner_count,tile_size);
    fprintf('done.\n');
    % return resulting features
    M=M(:,any(M)); % remove 0 elements
    Merkmale=M;
    
    %% show corners if specified
    if do_plot
        figure;
        title('Harris Detector using tiles and min distance');
        imshow(Img)
        hold on;            
        plot(M(1,:),M(2,:),'rx','MarkerSize',5);  
        hold off; 
    end
end

