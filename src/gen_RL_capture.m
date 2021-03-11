function gen_RL_capture(input_psf_path, file_name, sample_path)
%% generate capture in the target directory
%  resize the 3D volume
%  last update: 4/10/2020. YZ

% input_psf_path = '..\\..\\..\\psf';
% file_name = 'PSF_M_10_NA_0.5_d_750.0_from_-100_to_100_zspac_10_Nnum_7_OSR_3.mat';
% load PSF
load(sprintf('%s\\%s', input_psf_path, file_name));

% load sample in target path
% sample_path = '..\\trainning_dataset';
file_pattern = fullfile(sample_path, '*.tif');
tif_file = dir(file_pattern);

output_path1 = sprintf('%s\\LMF_recon', sample_path);
output_path2 = sprintf('%s\\LMF_cap', sample_path);
output_path3 = sprintf('%s\\high_res_volume_crop', sample_path);
mkdir(output_path1);
mkdir(output_path2);
mkdir(output_path3);
%% parameters
LFM_run = 3;
gpuFlag = true;
%% run
% for all files
for k = 1 : length(tif_file)

    base_filename = tif_file(k).name;
    stack = double(loadtiff(sprintf('%s\\%s',sample_path, base_filename)));
    stack = stack - min(stack(:));
    stack = stack / max(stack(:));

    
    % mod and resize
    stack = modcrop(stack, Nnum);
    stack = imresize3(stack, [size(stack, 1), size(stack, 2), size(H, 5)]);
    saveastiff(im2uint16(stack), sprintf('%s\\%s', output_path3, base_filename))
    
%     LFM capture
%     derectly run forward model
    xsize = [size(stack,1), size(stack,2)];
    msize = [size(H,1), size(H,2)];
    mmid = floor(msize/2);
    exsize = xsize + mmid;  % to make the size as 2^N after padding
    exsize = [ min( 2^ceil(log2(exsize(1))), 128*ceil(exsize(1)/128) ), min( 2^ceil(log2(exsize(2))), 128*ceil(exsize(2)/128) ) ];    
    zeroImageEx = gpuArray(zeros(exsize, 'single'));

    LFM_cap = forwardProjectGPU(H, stack, zeroImageEx, exsize);
    LFM_cap = gather(LFM_cap);
    LFM_cap = LFM_cap / max(LFM_cap(:));
    saveastiff(im2uint16(LFM_cap), sprintf('%s\\%s', output_path2, base_filename))

    % LFM recon
    Xguess = recon_LFM(LFM_cap, H, Ht, LFM_run, gpuFlag, []);  
    Xguess = gather(Xguess);
    saveastiff(im2uint16(Xguess/ max(Xguess(:))), sprintf('%s\\%s', output_path1, base_filename))
end
end

%% utility functions
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
end


function imgs = modcrop(imgs, modulo)
% imgs as a cell array, since different image might have different size.
    sz = size(imgs);
    sz = sz - mod(sz, modulo);
    imgs = imgs(1:sz(1), 1:sz(2), :);
end