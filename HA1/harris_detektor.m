function  Merkmale=harris_detektor(Image,varargin) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

 %% Input parser
    P = inputParser;

    % Liste der notwendigen Parameter
    % Ein Bild als Input ist zwingend notwendig
    P.addRequired('Image', @(x)isnumeric(x));

    % Liste der optionalen Parameter
    % Gr��e des Bildsegments
    % (default) = 10;
    P.addOptional('segment_length', 10, @(x) x ~=0);
    % Ecken- und Kantengewichtung mit k
    % (default) = 0.05;
    P.addOptional('k', 0.05, @(x)isnumeric(x));
    % Schwellenwert tau
    % (default) = 1;
    P.addOptional('tau', 1, @(x)isnumeric(x));
    P.addOptional('do_plot', 1, @(x)isnumeric(x));
    P.addOptional('tile_size', 4, @(x)isnumeric(x));
    P.addOptional('N', 20, @(x)isnumeric(x));
    P.addOptional('min_dist', 30, @(x)isnumeric(x));
    % Lese den Input
    P.parse(Image,varargin{:});

    %% Extrahiere die Variablen aus dem Input-Parser
    Img     = P.Results.Image; 
    % Img     = Bild;    % Auch moeglich, da die Variable Bild der Funtion direkt uebergeben wurde  
    slength   = P.Results.segment_length;
    k   = P.Results.k;
    tau = P.Results.tau;
    do_plot = P.Results.do_plot;
    sigma=1;
    %% Harris detector
    %% Convert to greyscale and perform sobel filter
    Img_gray=rgb_to_gray(Img);
    [Fx,Fy]=sobel_xy(Img_gray);
    %% Create Covariance matrix for each pixel to approximate Harris matrix
    G=zeros(size(Img,1),size(Img,2),4);
    H=zeros(size(Img,1),size(Img,2));
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
        G(:,:,i) = conv2(G(:,:,i),Weight_filter_1d,'same');
        G(:,:,i) = conv2(G(:,:,i),Weight_filter_1d','same');
    end
    % Calculate H Value for each pixel matrix, first determinant, then
    % trace
    H=G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3)-k*(G(:,:,1)+G(:,:,4)).^2;
    %H=(H>-tau) .* (H < tau); %Textureless
    %H=H<-tau; %Kante
    H=H>tau;
    %% show corners if specified
    if do_plot
        Img(:,:,1)=Img(:,:,1)+uint8(H*255);
        imshow(Img)
    end
  
  
  
  
  
  
end

