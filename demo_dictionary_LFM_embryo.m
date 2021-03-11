clc, clear
close all
%% this file would show how the learned dictionary would be used to 
%  learning.
%  last update: 4/10/2020. YZ

addpath('src')
addpath(genpath('ext'))
%% run scale up
Nnum = 13;
ksvd_conf.Tdata = 5; 

TV_lambda= 1;
TV_maxIter = 10;

% low res image
% dictionary_name = 'Learned_dictionary_size_1000_Tdata_5_peak_7';
% lores_stack = loadtiff('Drosophila_by_RL\\iter_2.tiff');
% output_path = sprintf('Drosophila_by_dictionary_size_1000_Tdata_5_peak_7_TV_lambda_%g', TV_lambda);


dictionary_name = 'data\\Learned_dictionary_size_1500_Tdata_5_peak_7_cut_0.01_overlap_4';
lores_stack = loadtiff('data\\test_bio_sample\\recon_40x_embryo_stack.tif');
lores_stack = double(lores_stack);
output_path = '.';
mkdir(output_path)
TV_enable = false;
%% scale up
for test_depth = 1: 6
    test_depth
    % load specific configuration for different depths
    conf = importdata(sprintf('%s\\depth_%d.mat', dictionary_name, test_depth ));
    
    % sample
    lores = lores_stack(:, :, test_depth);
    lores = double(lores);
    lores = lores / max(lores_stack(:));
    mean_lores = mean(lores(:));
    lores = modcrop({lores}, Nnum); % here since the low res image has a crop

    features = collect(conf, lores, conf.scale, conf.filters);
    features = double(features);

    % Encode features using OMP algorithm
    coeffs = omp(conf.dict_lores, conf.V_pca' * features, [], ksvd_conf.Tdata);
    patches = conf.dict_hires * coeffs; 
    % Reconstruct using patches' dictionary

    % Add low frequencies to each reconstructed patch
    patches = patches + collect(conf, lores, conf.scale, {});

    % Combine all patches into one image
    img_size = size(lores{1});
    grid = sampling_grid(img_size, conf.window, conf.overlap, conf.border, conf.scale);
    result = overlap_add(patches, img_size, grid);
    % TV desnoise
    
    % TV desnoise
    if TV_enable
        result = TV_2d_denoiser(result, TV_lambda, TV_maxIter);
    end
    result = result / mean(result(:)) * mean_lores;
    out_stack(:, :, test_depth) = result;
    
    saveastiff(im2uint16(result / max(result(:))), sprintf('%s\\depth_%d.tiff', output_path, test_depth))
    fprintf('.');
    figure(101), imshow(result, [])
end
saveastiff(im2uint16(out_stack), 'embryo.tif')
%% reconstruction
function imgs = modcrop(imgs, modulo)
% imgs as a cell array, since different image might have different size.
for i = 1:numel(imgs)
    sz = size(imgs{i});
    sz = sz - mod(sz, modulo);
    imgs{i} = imgs{i}(1:sz(1), 1:sz(2));
end
end