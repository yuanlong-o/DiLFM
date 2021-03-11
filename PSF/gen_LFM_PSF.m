function [H, Ht] = gen_LFM_PSF(pixelPitch, OSR, IMGSIZE_REF, M, Nnum, zmin, zspacing,...
    zmax, fobj, NA, lambda, d, n, fml)
%% light field PSF calculation
%  M: magnification
%  pixelPitch: pixel size of the sensor
%  OSR: optical super resolution ratio. For better accuracy?
%  Nnum: pixel number under each microlens
%  IMGSIZE_REF: half number of microlens used in this session
%  zmin, zmax, zspacing: axial sampling points distribution
%  fobj, focal length of objective
%  NA, numerical aperture
%  lambda, wavelength
%  d: focal length of ML     
%  n: refractive index
%  fml: focal length of microlens

%  last update: 11/8/2019. YZ

%% parameters
k = 2*pi*n/lambda; %% k
k0 = 2*pi*1/lambda; %% k

IMG_HALFWIDTH = max( Nnum*(IMGSIZE_REF + 1), 2*Nnum); % pixel number

% generate gloabal physical coordinate for all used microlens
x1space = (pixelPitch/OSR)*[-IMG_HALFWIDTH*OSR:1:IMG_HALFWIDTH*OSR]; % coordinate of total psf size
x2space = (pixelPitch/OSR)*[-IMG_HALFWIDTH*OSR:1:IMG_HALFWIDTH*OSR]; 
x1length = length(x1space);
x2length = length(x2space);

% generate local physical coordinate for a single microlens. 
x1MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2];
x2MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2];


% center point of image
centerPT = ceil(length(x1space)/2); 

% half of pixels in the sensor
halfWidth =  Nnum * OSR * (IMGSIZE_REF + 0 ); % half

x3objspace = [zmin:zspacing:zmax];
p3max = max(abs(x3objspace));
%% MLA phase calculation
% get the microlens phase, global one (stupid operation)
MLARRAY = calcML(fml, k0, x1MLspace, x2MLspace, x1space, x2space); 


%% PSF calculation
% sample space, physical coordinate under one pitch. Note without OSR
x1objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)]; 
x2objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)]; 
XREF = ceil(length(x1objspace)/2); % half length of the x1objspace
YREF = ceil(length(x1objspace)/2); % half length of the x1objspace

% get non-OSR coordinate
CP = ( (centerPT-1)/OSR+1 - halfWidth/OSR :1: (centerPT-1)/OSR+1 + halfWidth/OSR  );
H = zeros( length(CP), length(CP), length(x1objspace), length(x2objspace), length(x3objspace) );

% PSF matrix initialization
LFpsfWAVE_STACK = zeros(x1length, x2length, length(x3objspace)); % Note this is with OSR
psfWAVE_STACK = zeros(x1length, x2length, length(x3objspace));

% step 1, calculate the in-axis PSF

for eachpt=1 : length(x3objspace)      
    p1 = 0;
    p2 = 0;
    p3 = x3objspace(eachpt); % note the original code just needs non zeros one
    
    IMGSIZE_REF_IL = ceil(IMGSIZE_REF*( abs(p3)/p3max)); % linear?
    halfWidth_IL =  max(Nnum*OSR*(IMGSIZE_REF_IL + 0 ), 2*Nnum*OSR); % larger than one lens?
    % get a calculation area, 1d array
    centerArea_IL = max((centerPT - halfWidth_IL),1) :...
        min((centerPT + halfWidth_IL),length(x1space)); 
    
    disp(['size of center area = ' num2str(length(centerArea_IL)) 'X' num2str(length(centerArea_IL)) ]);    
    
    % call PSF calculator
    [psfWAVE LFpsfWAVE] = calcPSF(p1, p2, p3, fobj, NA, x1space, x2space, ...
        pixelPitch/OSR, lambda, MLARRAY, d, M, n,  centerArea_IL);
    
    psfWAVE_STACK(:,:,eachpt)  = psfWAVE;
    LFpsfWAVE_STACK(:,:,eachpt)= LFpsfWAVE;    

end 

% step 2, local shift the Debye PSF, and apply microlens array.
for i=1:length(x1objspace)*length(x2objspace)*length(x3objspace)
    % for all the point in space. note (a, b, c) is just index, have
    % nothing to do with real position. thus shift is validate, since we
    % shift the corresponding grid
    [a, b, c] = ind2sub([length(x1objspace) length(x2objspace) length(x3objspace)], i);  
    
    % read record Debye PSF, note this is full scale PSF
     psfREF = psfWAVE_STACK(:,:,c);  
     
     % shift based on pixel. if a is the center, no shift
     % otherwise, shift a 
     psfSHIFT= im_shift2(psfREF, OSR*(a-XREF), OSR*(b-YREF) );
     
     % apply fresnel propagation
     [f1,dx1,x1]=fresnel2D(psfSHIFT.*MLARRAY, pixelPitch/OSR, d,lambda);
     
     % shift back, this is quite interesting
     f1= im_shift2(f1, -OSR*(a-XREF), -OSR*(b-YREF) ); % shift back the pattern, why??
     % get image size
     xmin =  max( centerPT  - halfWidth, 1);
     xmax =  min( centerPT  + halfWidth, size(f1,1) );
     ymin =  max( centerPT  - halfWidth, 1);
     ymax =  min( centerPT  + halfWidth, size(f1,2) );
      
     % record, cut the edge
     f1_AP = zeros(size(f1));
     f1_AP( (xmin:xmax), (ymin:ymax) ) = f1( (xmin:xmax), (ymin:ymax) );
     % binding, since OSR
     [f1_AP_resize, x1shift, x2shift] = pixelBinning(abs(f1_AP.^2), OSR);           
     % slightly adjust the position
     f1_CP = f1_AP_resize( CP - x1shift, CP-x2shift ); 
     
     % record large PSF matrix
     H(:,:,a,b,c) = f1_CP; % here H is a just the mapping function under one pitch, contain everything
     
end

% Ht calculation
Ht = calcHt(H);
%% utility functions
function MLARRAY = calcML(fml, k, x1MLspace, x2MLspace, x1space, x2space)
%% microlens array phase calculation
% input:
%  fml: focal length of microlens
%  k: wave vector
%  x1MLspace: local physical coordinate for a single microlens, OSR within
%  x1space: globale physical coordinate, OSR within
% output:
%

% function for generating microlens array
x1length = length(x1space);
x2length = length(x2space);
x1MLdist = length(x1MLspace); % OSR with single lenslet pixels
x2MLdist = length(x2MLspace);
x1center = find(x1space==0); % stupid code, if not?
x2center = find(x2space==0);

% get the microlens center
x1centerALL = [  (x1center: -x1MLdist:1)  (x1center + x1MLdist: x1MLdist :x1length)];
x1centerALL = sort(x1centerALL);
x2centerALL = [  (x2center: -x2MLdist:1)  (x2center + x2MLdist: x2MLdist :x2length)];
x2centerALL = sort(x2centerALL);

% local phase
patternML = zeros( length(x1MLspace), length(x2MLspace) );

% loop for each pixel under a microlens
for a=1:length(x1MLspace),
    for b=1:length(x2MLspace),        
        x1 = x1MLspace(a);
        x2 = x2MLspace(b);
        xL2norm = x1^2 + x2^2;
        
        patternML(a,b) = exp(-i*k/(2*fml)*xL2norm);   
    end
end

% global
MLspace = zeros( length(x1space), length(x2space) );
MLcenters = MLspace;
for a=1:length(x1centerALL),
    for b=1:length(x2centerALL),
        MLcenters( x1centerALL(a), x2centerALL(b)) = 1; % global center
    end
end
% convolution to allot the phase
MLARRAY = conv2(MLcenters, patternML, 'same');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function new_im = im_shift2(img, SHIFTX, SHIFTY)
%% shift the image based on shift vector SHIFTX and SHIFTY
% calculate the image shift, without consideration of boundary
eqtol = 1e-10;

xlength = size(img,1);
ylength = size(img,2);

if abs(mod(SHIFTX,1))>eqtol | abs(mod(SHIFTY,1))>eqtol
   error('SHIFTX and SHIFTY should be integer numbers');
end

% if SHIFTX >= xlength | SHIFTY >= ylength,
%    error('SHIFTX  and SHIFTY should be smaller than size(img,1) and size(img,2), respectively');
% end

SHIFTX = round(SHIFTX);
SHIFTY = round(SHIFTY);

new_im = zeros(xlength, ylength, size(img,3) );

if SHIFTX >=0 & SHIFTY >= 0,
    new_im( (1+SHIFTX:end), (1+SHIFTY:end),:) = img( (1:end-SHIFTX), (1:end-SHIFTY),:);
elseif SHIFTX >=0 & SHIFTY < 0,
    new_im( (1+SHIFTX:end), (1:end+SHIFTY),:) = img( (1:end-SHIFTX), (-SHIFTY+1:end),:);
elseif SHIFTX <0 & SHIFTY >= 0,
    new_im( (1:end+SHIFTX), (1+SHIFTY:end),:) = img( (-SHIFTX+1:end), (1:end-SHIFTY),:);
else
    new_im( (1:end+SHIFTX), (1:end+SHIFTY),:) = img( (-SHIFTX+1:end), (-SHIFTY+1:end),:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [psf, LFpsf] = calcPSF(p1, p2, p3, fobj, NA, x1space, x2space, scale, lambda, MLARRAY, fml, M, n, centerArea)
%% Debye PSF calculator
%  centerarea: the response area under such depth. OSR within

%  output:
%  psf: the output PSF from Debye model
%  LFpsf: corresponding LF psf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
k = 2*pi*n/lambda;
alpha = asin(NA/n);
x1length = length(x1space);
x2length = length(x2space);
zeroline = zeros(1, length(x2space) );
 
% global PSF
pattern = zeros(x1length, x2length);

% center pixel
centerPT = ceil(length(x1space)/2);

% step 1, Debye PSF
for a = centerArea(1) : centerPT % 1d. a and b is just index
    patternLine = zeroline; % initialization.
    for b = a:centerPT % 2d
        % we could see the program would only update a small area
        x1 = x1space(a);
        x2 = x2space(b);          
        xL2normsq = (((x1+M*p1)^2+(x2+M*p2)^2)^0.5)/M;        
       
        v = k*xL2normsq*sin(alpha);    
        u = 4*k*(p3*1)*(sin(alpha/2)^2);
        Koi = M/((fobj*lambda)^2)*exp(-1i*u/(4*(sin(alpha/2)^2)));
        intgrand = @(theta) (sqrt(cos(theta))) .* (1+cos(theta))  .*  (exp(-(1i*u/2)* (sin(theta/2).^2) / (sin(alpha/2)^2)))  .*  (besselj(0, sin(theta)/sin(alpha)*v))  .*  (sin(theta));
        I0 = integral(@(theta)intgrand (theta),0,alpha);  
                
        patternLine(1,b) =  Koi*I0; % interesting
    end
    pattern(a,:) = patternLine;    % this is just 1/8 of the whole PSF
end

% rotate and copy
patternA = pattern( (1:centerPT), (1:centerPT) );
patternAt = fliplr(patternA);

pattern3D = zeros(size(pattern,1), size(pattern,2), 4);
pattern3D(:,:,1) = pattern;
pattern3D( (1:centerPT), (centerPT:end),1 ) = patternAt;
pattern3D(:,:,2) = rot90( pattern3D(:,:,1) , -1);
pattern3D(:,:,3) = rot90( pattern3D(:,:,1) , -2);
pattern3D(:,:,4) = rot90( pattern3D(:,:,1) , -3);
pattern = max(pattern3D,[],3); % maximum projection

% step 2, response in sensor plane
[f1,dx1,x1]=fresnel2D(pattern.*MLARRAY,scale,1*fml,lambda);
psf = pattern; % note this is field, not intensity
LFpsf = f1;

%% Fresnel propagation
function [f1,dx1,x1]=fresnel2D(f0,dx0,z,lambda)
Nx = size(f0,1);
Ny = size(f0,2);
k = 2*pi/lambda;

%
du = 1./(Nx*dx0);
u = [0:ceil(Nx/2)-1 ceil(-Nx/2):-1]*du; 
dv = 1./(Ny*dx0);
v = [0:ceil(Ny/2)-1 ceil(-Ny/2):-1]*dv; 

H = exp(-i*2*pi^2*(repmat(u',1,length(v)).^2+repmat(v,length(u),1).^2)*z/k); 
f1 = exp(i*k*z)*ifft2( fft2(f0) .* H ); 
dx1 = dx0;
x1 = [-Nx/2:Nx/2-1]*dx1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pixel binding to cancel the OSR
function [OIMG, x1shift, x2shift] = pixelBinning(SIMG, OSR)
% Input:
% SIMG: super resolved input image

% do imaging pixel binning
x1length = size(SIMG,1);
x2length = size(SIMG,2);

x1center = (x1length-1)/2 + 1;
x2center = (x2length-1)/2 + 1;

x1centerinit = x1center - (OSR-1)/2;
x2centerinit = x2center - (OSR-1)/2;
x1init = x1centerinit -  floor(x1centerinit/OSR)*OSR ;
x2init = x2centerinit -  floor(x2centerinit/OSR)*OSR ;


x1shift = 0;
x2shift = 0;
if x1init<1,
    x1init = x1init + OSR;
    x1shift = 1;
end
if x2init<1,
    x2init = x2init + OSR;
    x2shift = 1;
end
halfWidth = length( (x1init:x1center-1) );
SIMG_crop = SIMG( [ (x1init:x1center-1) x1center x1center+1:x1center+halfWidth ],  [ (x2init:x2center-1) x2center x2center+1:x2center+halfWidth ] );

[m,n]=size(SIMG_crop); %M is the original matrix
SIMG_crop = sum( reshape(SIMG_crop,OSR,[]) ,1 );
SIMG_crop=reshape(SIMG_crop,m/OSR,[]).'; %Note transpose
SIMG_crop=sum( reshape(SIMG_crop,OSR,[]) ,1);
OIMG =reshape(SIMG_crop,n/OSR,[]).'; %Note transpose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Ht = calcHt(H)
%% calculate backward projection matrix Ht
% while Ht has the same dimension of H
% Note that H is not fully represented, thus they take a interesting
% calculation strategy to calculate Ht


%  depth by depth
Hsize = size(H);
Hsize1 = Hsize(1);
Nnum = Hsize(3); % number of virtural pixel number of a microlens
x3length = Hsize(5);

tmpsize = ceil(size(H,1)/Nnum); % calculate roughly how many microlens in this psf calibration
if mod(tmpsize,2) == 1,
    imgsize = (tmpsize+2)*Nnum;
else
    imgsize = (tmpsize+3)*Nnum, % calculte image size in terms of pixel number 
end

zeroprojection = zeros(imgsize, imgsize);
imcenter = ceil(imgsize/2); % center of the image
imcenterinit = imcenter - ceil(Nnum/2); 

Ht = zeros(Hsize);
for aa=1:1:Nnum, % run over the virtral pixels
    for bb=1:1:Nnum,
        temp = zeroprojection;
        temp( imcenterinit+aa, imcenterinit+bb ) = 1; % give 1 to the corresponding pixel

        tempback = backwardProject(H, temp, Nnum) ; % do backward project, now tempback would be the weight matrix in sample side
        tempback_cut = tempback( ( imcenter - (Hsize1-1)/2 - 0*Nnum :imcenter + (Hsize1-1)/2 + 0*Nnum), ( imcenter - (Hsize1-1)/2 - 0*Nnum :imcenter + (Hsize1-1)/2 + 0*Nnum) , :);
        % simply cut to H size
        for cc=1:x3length,
            tempback_shift(:,:,cc) = im_shift2(tempback_cut(:,:,cc), (ceil(Nnum/2)-aa) , (ceil(Nnum/2)-bb) ); % why do shift here?
        end
        % shift because later would conduct convolution, so make sure the
        % PSF in the center
        Ht(:,:,aa,bb,:) = tempback_shift; % well this program is confusing, (aa, bb) should be object space coordinate, here they put it as sensor space
    end       
end


function Backprojection = backwardProject(H, projection, Nnum )
% run back propgation
% Input:
% H: 5 dim matrix, PSF. First 2 dims are pixels over the WHOLE SENSOR, last three are
% spatial voxels restricted by a single microlens
% projection: 2d matrix, indicate which pixels in the sensor plane would
% like to be back propogated
% Nnum: virtral pixel number corresponding to a single microlens
% Output
% Backprojection: 3d matrix, give the 

% Commented by YZ

x3length = size(H,5); % voxel depth dimension
Backprojection = zeros(size(projection, 1), size(projection, 1), x3length );


for aa=1:Nnum,
    for bb=1:Nnum,
        for cc=1:x3length, % run over the pixels in 3D space
                      
            Ht = imrotate( H(:,:,aa,bb,cc), 180); % rotate the psf generated by the specific voxel
            % equal to H(end : -1 : 1, end : -1 : 1, aa, bb, cc)
            tempSlice = conv2(projection, Ht, 'same'); % 2D convolution, between this PSF and sensor active pattern?
            Backprojection((aa:Nnum:end) , (bb:Nnum:end),cc) = ...
                Backprojection((aa:Nnum:end) , (bb:Nnum:end),cc) + tempSlice( (aa:Nnum:end) , (bb:Nnum:end) ); % this part is confusing, and I think this is wrong. this should not be aa, bb but u and v
            % periodical summation. Note the result is accumulate
            % briliant method
        end
    end
end

