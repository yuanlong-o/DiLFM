%% here recorld all the functions that needs for phase space LFM reconstruction
%  this one is exactly what the theory suggested

%  Last update: 11/16/2019
function Xguess = recon_theory_phase_space_LFM(LFIMG, H, Ht, weight, maxIter, gpuFlag)
Nnum = size(H, 3);
%% define operator
if gpuFlag == true
    xsize = [size(LFIMG,1), size(LFIMG,2)];
    msize = [size(H,1), size(H,2)];
    mmid = floor(msize/2);
    exsize = xsize + mmid;  % to make the size as 2^N after padding
    exsize = [ min( 2^ceil(log2(exsize(1))), 128*ceil(exsize(1)/128) ), min( 2^ceil(log2(exsize(2))), 128*ceil(exsize(2)/128) ) ];    
    zeroImageEx = gpuArray(zeros(exsize, 'single'));
    % push to GPU
%     Ht = gpuArray(Ht);
%     H = gpuArray(H);
%     LFIMG = gpuArray(LFIMG);
    
    backwardFUN = @(psft, projection, u, v, Nnum) backwardProjectGPU_phasespace_theory(psft, gpuArray(projection), zeroImageEx, exsize, u, v, Nnum);
    forwardFUN = @(psf, Xguess, Nnum) forwardProjectGPU_phasespace_theory(psf, gpuArray(Xguess), zeroImageEx, exsize, Nnum); % one use H and one use Ht
else
	forwardFUN =  @(psf, Xguess, Nnum) forwardProjectACC_phasespace_theory(H, Xguessv, Nnum); % build the function: forward and backward
    backwardFUN = @(psft, projection, u, v, Nnum) backwardProject_phasespace_theory(psft,  projection, u, v, Nnum); 
end


% Xguess = contrastAdjust(Xguess, contrast);

% update index
[index1, index2] = gen_spiral(Nnum);

% initialization
Xguess=ones(size(LFIMG,1),size(LFIMG,2),size(H,5));
Xguess=gpuArray(Xguess./sum(Xguess(:)).*sum(LFIMG(:))); % this is a poor initialize
% generate Htf
% Htf = zeros(size(LFIMG,1),size(LFIMG,2), size(H,5), Nnum ,Nnum);
% for u=1:Nnum
%     for v=1:Nnum
%         if weight(u,v)==0
%             continue;
%         else
%             Htf(:,:,:,u,v)= backwardFUN(squeeze(Ht(:,:,u,v,:)), gpuArray(ones(size(LFIMG))));
%         end
%     end
% end


% main iteration
for i=1:maxIter
    tic;
    for u_2=1:Nnum
        for v_2=1:Nnum 
            u=index1((u_2-1)*Nnum+v_2); % marginal to center, rebundant information first
            v=index2((u_2-1)*Nnum+v_2);
%             u = u_2;
%             v = v_2;
            if weight(u,v)==0 % update weight
                continue;
            else
                HXguess = forwardFUN(gpuArray(squeeze(H(:,:,u,v,:))), Xguess, Nnum);  % Note HXguess would be small
%                 figure, imagesc(HXguess)
                % selective
                HXguess = squeeze(LFIMG(u : Nnum : end, v : Nnum : end))./HXguess;
                HXguess(~isfinite(HXguess))=0;           
                
                XguessCor = backwardFUN(gpuArray(squeeze(Ht(:,:,u,v,:))),HXguess, u, v, Nnum) ;  
                buf = backwardFUN(gpuArray(squeeze(Ht(:,:,u,v,:))), gpuArray(ones(size(LFIMG)/Nnum)), u, v, Nnum); % without pre-computing, why this?
                XguessCor=Xguess.*XguessCor./buf;
%                 XguessCor=Xguess.*XguessCor;
                XguessCor(find(isnan(XguessCor))) = 0;
                XguessCor(find(isinf(XguessCor))) = 0;
                XguessCor(XguessCor<0 ) = 0;
                % update
                Xguess=XguessCor.*weight(u,v)+(1-weight(u,v)).*Xguess;
%                 Xguess=XguessCor.*weight(u,v)+ Xguess;
                Xguess(Xguess<0) = 0;              
            end
        end
    end
    ttime1=toc;
    disp(['iter ',num2str(i),' | ',num2str(maxIter),', phase-space deconvolution took ',num2str(ttime1),' secs']);
end

%% utility function
%% generate spiral index
function [i_index, j_index] = gen_spiral(Nnum)

if mode(Nnum) == 0
    loop_n = Nnum / 2;
else
    loop_n = (Nnum + 1) / 2;
end
i_index = zeros(1, Nnum * Nnum);
j_index = zeros(1, Nnum * Nnum);

start = 1;
for k = 1 :loop_n
    current_size = Nnum - 2 * (k - 1);
    bias = k - 1;
    if current_size > 1
        i_new = [1 : current_size , current_size * ones(1, current_size - 1), ...
            current_size-1 : -1 : 1, ones(1, current_size - 2)];
        
        j_new = [ones(1, current_size), 2 : current_size, ...
            current_size * ones(1, current_size - 1), current_size - 1 : -1 : 2];
        
        i_index(start : start + 4 * (current_size - 1) - 1) = i_new + bias;
        j_index(start : start + 4 * (current_size - 1) - 1) = j_new + bias;
        start = start + 4 * (current_size - 1);
    else
        i_new = 1;
        j_new = 1;
        i_index(start) = i_new + bias;
        j_index(start) = j_new + bias;        
    end
   
end

%% forward projection function
function projection_out = forwardProjectACC_phasespace_theory(H, realspace, Nnum)
% build forward function
projection=zeros(size(realspace,1),size(realspace,2));
for z=1:size(realspace,3)
    projection=projection+conv2(realspace(:,:,z),H(:,:,z),'same');
end
projection_out = projection((Nnum + 1)/2 : Nnum : end, (Nnum + 1)/2 : Nnum : end);

%% backprojection function
function Backprojection = backwardProject_phasespace_theory(Ht,projection, u, v, Nnum)
big_projection = zeros(Nnum * size(projection));
big_projection(u : Nnum : end, v : Nnum : end) = projection;
Backprojection=zeros(size(big_projection,1),size(big_projection,2),size(Ht,3));
for z=1:size(Ht,3)
    Backprojection(:,:,z)=conv2(big_projection,Ht(:,:,z),'same');
end

%% backprojection function with GPU
function Backprojection = backwardProjectGPU_phasespace_theory(Ht, projection, zeroImageEx, exsize, u, v, Nnum)
% generate backward projection 
big_projection = gpuArray.zeros(Nnum * size(projection));
big_projection(u : Nnum : end, v : Nnum : end) = projection;
Backprojection = gpuArray.zeros(size(big_projection,1),size(big_projection,2),size(Ht,3));

for z=1:size(Ht,3)
    Backprojection(:,:,z)=conv2FFT(big_projection, Ht(:,:,z), zeroImageEx, exsize);
end

%% forward function with GPU
function projection_out = forwardProjectGPU_phasespace_theory( H, realspace, zeroImageEx, exsize, Nnum)
projection=gpuArray.zeros(size(realspace,1),size(realspace,2));
for z=1:size(realspace,3)
    projection=projection+conv2FFT(realspace(:,:,z),H(:,:,z), zeroImageEx, exsize);
end
projection_out = projection((Nnum + 1)/2 : Nnum : end, (Nnum + 1)/2 : Nnum : end); % additional selection
