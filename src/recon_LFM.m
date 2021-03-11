%% here recorld all the functions that needs for LFM reconstruction
%  Last update: 11/15/2019
function Xguess = recon_LFM(LFIMG, H, Ht, maxIter, gpuFlag, output_path)
%% define operator
if gpuFlag == true
    xsize = [size(LFIMG,1), size(LFIMG,2)];
    msize = [size(H,1), size(H,2)];
    mmid = floor(msize/2);
    exsize = xsize + mmid;  % to make the size as 2^N after padding
    exsize = [ min( 2^ceil(log2(exsize(1))), 128*ceil(exsize(1)/128) ), min( 2^ceil(log2(exsize(2))), 128*ceil(exsize(2)/128) ) ];    
    zeroImageEx = gpuArray(zeros(exsize, 'single'));
    
    backwardFUN = @(projection) backwardProjectGPU(Ht, projection, zeroImageEx, exsize);
    forwardFUN = @(Xguess) forwardProjectGPU( H, Xguess, zeroImageEx, exsize); % one use H and one use Ht
else
	forwardFUN =  @(Xguess) forwardProjectACC( H, Xguess, CAindex ); % build the function: forward and backward
    backwardFUN = @(projection) backwardProjectACC(Ht, projection, CAindex ); 
end


% Determined initial guess
Htf = backwardFUN(LFIMG); 
Xguess = Htf;
% Xguess = contrastAdjust(Xguess, contrast);
for i=1:maxIter
    tic;
    HXguess = forwardFUN(Xguess);
    HXguessBack = backwardFUN(HXguess);
    errorBack = Htf./HXguessBack;
    Xguess = Xguess.*errorBack; 
    Xguess(find(isnan(Xguess))) = 0;
    ttime = toc;
    disp(['  iter ' num2str(i) ' | ' num2str(maxIter) ', took ' num2str(ttime) ' secs']);
    
    Xguess = gather(Xguess);
%     save(sprintf('recon_%s_iter_%d_lambda.mat', LFM_img_name, maxIter), 'Xguess')
    saveastiff(im2uint16(Xguess / max(Xguess(:))), sprintf('%s\\recon_iter_%d.tiff', output_path, i))
end
%% utility function
%% forward projection function
function TOTALprojection = forwardProjectACC( H, realspace, CAindex)
% build forward function
Nnum = size(H,3);
zerospace = zeros(  size(realspace,1),   size(realspace,2), 'single');
TOTALprojection = zerospace;
% just loop for a small space behind a microlens
for aa=1:Nnum
    for bb=1:Nnum
        for cc=1:size(realspace,3) % for each spatial points

            % read out the PSF
            Hs = squeeze(H( CAindex(cc,1):CAindex(cc,2), CAindex(cc,1):CAindex(cc,2) ,aa,bb,cc));          
            tempspace = zerospace;
            tempspace( (aa:Nnum:end), (bb:Nnum:end) ) = realspace( (aa:Nnum:end), (bb:Nnum:end), cc); % get the pixels from the same place out
            projection = conv2(tempspace, Hs, 'same');
            TOTALprojection = TOTALprojection + projection;            
        end
    end
end


%% backprojection function
function Backprojection = backwardProjectACC(Ht, projection, CAindex )
x3length = size(Ht,5);
Nnum = size(Ht,3);
Backprojection = zeros(size(projection, 1), size(projection, 2), x3length);
zeroSlice = zeros(  size(projection,1) , size(projection, 2));

for cc=1:x3length
    tempSliceBack = zeroSlice;
    for aa=1:Nnum
        for bb=1:Nnum
                          
            Hts = squeeze(Ht( CAindex(cc,1):CAindex(cc,2), CAindex(cc,1):CAindex(cc,2) ,aa,bb,cc));                  
            tempSlice = zeroSlice;
            tempSlice( (aa:Nnum:end) , (bb:Nnum:end) ) = projection( (aa:Nnum:end) , (bb:Nnum:end) );
            tempSliceBack = tempSliceBack + conv2(tempSlice, Hts, 'same');   

        end
    end
    Backprojection(:,:,cc) = Backprojection(:,:,cc) + tempSliceBack;
end

function Xguess_out  = contrastAdjust(Xguess_in, Contrast)

XguessMax = max(Xguess_in(:));
Xguess = Xguess_in/XguessMax;
Xguess = Contrast*(2*Xguess - 1)+1;
Xguess_out = 0.5*XguessMax*Xguess;

%% backprojection function with GPU
function Backprojection = backwardProjectGPU(Ht, projection, zeroImageEx, exsize)
% generate backward projection 

Nnum = size(Ht,3);
x3length = size(Ht,5);
Backprojection = gpuArray.zeros(size(projection, 1), size(projection, 2), x3length , 'single');
zeroSlice = gpuArray.zeros(size(projection,1) , size(projection, 2) , 'single');


for cc=1:x3length,
    tempSliceBack = zeroSlice;
    for aa=1:Nnum,
        for bb=1:Nnum,                  
            Hts = gpuArray(squeeze(Ht(:,:, aa,bb,cc)));        
            tempSlice = zeroSlice;
            tempSlice( (aa:Nnum:end) , (bb:Nnum:end) ) = projection( (aa:Nnum:end) , (bb:Nnum:end) );
            tempSliceBack = tempSliceBack + conv2FFT(tempSlice, Hts, zeroImageEx, exsize);
        end
    end
    Backprojection(:,:,cc) = Backprojection(:,:,cc) + tempSliceBack;
end

%% forward function with GPU
function TOTALprojection = forwardProjectGPU( H, realspace, zeroImageEx, exsize)

Nnum = size(H,3);
zerospace = gpuArray.zeros(  size(realspace,1),   size(realspace,2), 'single');
TOTALprojection = zerospace;

for aa=1:Nnum
    for bb=1:Nnum
        for cc=1:size(realspace,3),
    
            Hs = gpuArray(squeeze(H( :,:,aa,bb,cc)));    
            tempspace = zerospace;
            tempspace( (aa:Nnum:end), (bb:Nnum:end) ) = realspace( (aa:Nnum:end), (bb:Nnum:end), cc);
            projection = conv2FFT(tempspace, Hs, zeroImageEx, exsize);
            TOTALprojection = TOTALprojection + projection;            

        end
    end
end
TOTALprojection = double(TOTALprojection);


