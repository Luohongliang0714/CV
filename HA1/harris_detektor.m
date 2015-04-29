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
    P.addOptional('tau', 400, @(x)isscalar(x)&isnumeric(x));
    P.addOptional('do_plot', 1, @(x)islogical(x));
    P.addOptional('tile_size', 200, @(x) (isscalar(x) | prod(size(x)==[2,1]))&isnumeric(x));
    P.addOptional('N', 2, @(x)isscalar(x)&isnumeric(x)&x>1);
    P.addOptional('min_dist', 3, @(x)isscalar(x)&isnumeric(x));
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
    if isscalar(tile_size)
        tile_size=[tile_size;tile_size];
    end
    N=P.Results.N;
    min_dist=P.Results.min_dist;
    sigma=1; % set to fixed value
    
    %% Harris detector
    %% Convert to greyscale and perform sobel filter
    Img_gray=rgb_to_gray(Img);
    [Fx,Fy]=sobel_xy(Img_gray);
    %% Create Jacobian matrix for each pixel to approximate Harris matrix
    G=zeros(size(Img,1),size(Img,2),4);
    G(:,:,1)=Fx.^2;
    G(:,:,2)=Fx.*Fy;
    G(:,:,3)=G(:,:,2);
    G(:,:,4)=Fy.^2;
    segment_half = floor(slength/2);            % Halbe Fesntergroesse
    ind = -segment_half:segment_half;               % indices (Bsp.: -1,0,1)
    gau = exp(-(ind.^2)/(2*sigma^2))';    % 1D Gauss
    Weight_filter = gau;                           % 2 Separable Filter ergeben das 2D Filter
    Weight_filter_1d = Weight_filter./sum(Weight_filter);                 % Summe aller Koeffizienten ergibt 1
    G(:,:,1) = conv2(double(G(:,:,1)),Weight_filter_1d,'same');
    G(:,:,1) = conv2(double(G(:,:,1)),Weight_filter_1d','same');
    G(:,:,2) = conv2(double(G(:,:,2)),Weight_filter_1d,'same');
    G(:,:,2) = conv2(double(G(:,:,2)),Weight_filter_1d','same');
    G(:,:,3) = G(:,:,2);
    G(:,:,4) = conv2(double(G(:,:,4)),Weight_filter_1d,'same');
    G(:,:,4) = conv2(double(G(:,:,4)),Weight_filter_1d','same');
    
    % Calculate H Value for each pixel matrix, first determinant, then
    % trace
    H=G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3)-k*(G(:,:,1)+G(:,:,4)).^2;
    %H=(H>-tau) .* (H < tau); %Textureless
    %H=H<-tau; %Kante
    H_mask=H>tau; % Ecke

    % free memory, G not needed anymore
    clearvars G;
    %% Purge all smaller features within min_dist
    % Delate image with [min_dist*2+1,min_dist*2+1] mask
    % --> search for max element within window
    % only allow unchanged values -> smaller ones are purged
    % purge features at the border within min_dist first
    H_border=padarray(ones(size(H,1)-min_dist*2,size(H,2)-min_dist*2),[min_dist,min_dist]);
    H=H.*H_border;
    % Shift image (2*min_dist+1)^2-1 times
    shift_mat=zeros(size(H,1),size(H,2),2);
    shift_mat(:,:,1)=H;
    for i=-min_dist:min_dist
        for j=-min_dist:min_dist
            if (j~=0 || i~=0)
                % shift matrix (borders dont care because they are 0 anyways)
                shift_mat(:,:,2)=circshift(H,[i,j]); % faster
                shift_mat(:,:,1)=max(shift_mat,[],3);
                % if the values are the same, delete one of the features,
                % quite slow
                H=H.*(shift_mat(:,:,2)~=H);
            end
        end
    end
    % Maximum where value did not change
    H_dil=(H==shift_mat(:,:,1));
    % Apply the filters
    H=H.*H_mask.*H_dil;
    
    %% Segment the image in tiles and take the N best features
    % Check if the tilesize is well chosen, if not only full tiles are used
    randy=mod(size(H,1),tile_size(2));
    randx=mod(size(H,2),tile_size(1));
    % number of fill tiles in x and y direction
    num_tilesy=((size(H,1)-randy)/tile_size(2));
    num_tilesx=((size(H,2)-randx)/tile_size(1));
    M=zeros(2,N*num_tilesy*num_tilesx); % preallocate for speed
    corner_count=0;
    for i=1:num_tilesy
        for j=1:num_tilesx
            % get the tile
            tile=H(1+(i-1)*tile_size(2):1+i*tile_size(2)-1,1+(j-1)*tile_size(1):1+j*tile_size(1)-1);
            % sort its entries
            [sortedValues,sortIndex] = sort(tile(:),'descend');
            if size(sortIndex,1)>=N % more entries than N found
                maxIndex = sortIndex(1:N);
                maxIndex(sortedValues(1:N)==0)=[]; % get the N largest, get rid of 0 elements
            else % less than N found
                maxIndex = sortIndex(1:end);
                maxIndex(sortedValues(1:end)==0)=[];
            end
            [row,col] = ind2sub(size(tile), maxIndex);
            row=row+(i-1)*tile_size(2);
            col=col+(j-1)*tile_size(1);
            if size(maxIndex,1)~=0 
                % add features to result matrix
                M(:,corner_count+1:corner_count+size(row,1))=[col,row]';
                corner_count=corner_count+size(row,1);
            end
        end
    end
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

