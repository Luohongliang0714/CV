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

    % free memory
    clearvars G;
    %reshape H matrix with tile_size
    if isscalar(tile_size)
        randy=mod(size(H,1),tile_size);
        randx=mod(size(H,2),tile_size);
        num_tiles=((size(H,1)-randy)/tile_size)*((size(H,2)-randx)/tile_size);
        %Tiles=reshape(H(1:size(H,1)-randy,1:size(H,2)-randx),tile_size,tile_size,num_tiles);
    else
        % falls vector
       % Tiles=reshape(H,tile_size(1),tile_size(2));
    end
    % Filter all Features which do not have the min_dist. 
    % Search for next_neighbor and keep largest one
    % filter image and check which features are close to others
%     [feat_r,feat_c]=find(H);
%     feat_mat=[feat_r,feat_c]';
%     dist=zeros(size(feat_mat));
%     for i=1:size(feat_r,1)
%         dist(1,:)=(feat_mat(1,:)-feat_mat(1,i)).^2;
%         dist(2,:)=(feat_mat(2,:)-feat_mat(2,i)).^2;
%         dist2=dist(1,:)+dist(2,:); % Wurzel brauchen wir nicht
%         dist2(dist2==0)=[];
%         [nn_val,nn_ind]=min(dist2);
%         if H((feat_mat(:,i)>H(feat_mat(nn_ind)))&nn_val<=min_dist)
%             H(feat_mat(nn_ind))=0;
%         end
%         fprintf('%3d%% \n\n',i/size(feat_r,1)*100);
%     end
    % Purge all smaller features within min_dist
    % Delate image with [min_dist*2+1,min_dist*2+1] mask
    % --> search for max element within window
    % only allow unchanged values -> smaller ones are purged
    % purge features at the border within min_dist
    H_border=padarray(ones(size(H,1)-min_dist*2,size(H,2)-min_dist*2),[min_dist,min_dist]);
    H=H.*H_border;
    shift_mat=zeros(size(H,1),size(H,2),2);
    shift_mat(:,:,1)=H;
    for i=-min_dist:min_dist
        for j=-min_dist:min_dist
            shift_mat(:,:,2)=circshift(H,[i,j]);
            shift_mat(:,:,1)=max(shift_mat,[],3);
        end
    end
    H_dil=(H==shift_mat(:,:,1));
    H=H.*H_mask.*H_dil;
    num_tilesy=((size(H,1)-randy)/tile_size);
    num_tilesx=((size(H,2)-randx)/tile_size);
    M=zeros(2,N*num_tilesy*num_tilesx);
    corner_count=0;
    for i=1:num_tilesy
        for j=1:num_tilesx
            tile=H(1+(i-1)*tile_size:1+i*tile_size-1,1+(j-1)*tile_size:1+j*tile_size-1);
            [sortedValues,sortIndex] = sort(tile(:),'descend');
            if size(sortIndex,1)>=N
                maxIndex = sortIndex(1:N);
                maxIndex(sortedValues(1:N)==0)=[];
            else 
                maxIndex = sortIndex(1:end);
                maxIndex(sortedValues(1:end)==0)=[];
            end
            [row,col] = ind2sub(size(tile), maxIndex);
            row=row+(i-1)*tile_size;
            col=col+(j-1)*tile_size;
            %b=[col,row]'
            if size(maxIndex,1)~=0 
                M(:,corner_count+1:corner_count+size(row,1))=[col,row]';
                corner_count=corner_count+size(row,1);
            end
        end
    end
    M=M(:,any(M));
    Merkmale=M;
%     feat_mat=M;
%     % remove elements too close
%     dist=zeros(size(feat_mat));
%     for i=1:size(M,2)
%         dist(1,:)=(feat_mat(1,:)-feat_mat(1,i)).^2;
%         dist(2,:)=(feat_mat(2,:)-feat_mat(2,i)).^2;
%         dist2=dist(1,:)+dist(2,:); % Wurzel brauchen wir nicht
%         dist2(dist2==0)=800000;
%         [nn_val,nn_ind]=min(dist2);
%         if (H(feat_mat(2,i),feat_mat(1,i))>H(feat_mat(2,nn_ind),feat_mat(1,nn_ind)))&&nn_val<=sqrt(min_dist)
%             fprintf('%d',H(feat_mat(2,nn_ind),feat_mat(1,nn_ind)));
%             H(feat_mat(2,nn_ind),feat_mat(1,nn_ind))=0;
%         end
%         fprintf('%3d%% \n\n',i/size(M,2)*100);
%     end
%     corner_count=0;
%     M_new=zeros(size(M));
%     for i=1:num_tilesy
%         for j=1:num_tilesx
%             tile=H(1+(i-1)*tile_size:1+i*tile_size,1+(j-1)*tile_size:1+j*tile_size);
%             [sortedValues,sortIndex] = sort(tile(:),'descend');
%             maxIndex = sortIndex(1:N);
%             maxIndex(sortedValues(1:N)==0)=[];
%             [row,col] = ind2sub(size(tile), maxIndex);
%             row=row+(i-1)*tile_size;
%             col=col+(j-1)*tile_size;
%             M_new(:,corner_count+1:corner_count+size(row,1))=[col,row]';
%             corner_count=corner_count+size(row,1);
%         end
%     end
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

