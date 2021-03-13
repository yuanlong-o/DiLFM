function out_stack = lfm_dictionary_test(input_psf_path, file_name, TV_config, dictionary_name, lores_stack, Tdata)
% for dictionary test
%  last udpate: 3/11/2021.
load(sprintf('%s/%s', input_psf_path, file_name));

% dictionary_name = 'Learned_dictionary_size_1500_Tdata_5_peak_4_cut_0.01_overlap_9';
assert(size(lores_stack, 3) == size(H, 5))
% mkdir(output_path)
%% scale up
for test_depth = 1: 21
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
    coeffs = omp(conf.dict_lores, conf.V_pca' * features, [], Tdata);
    patches = conf.dict_hires * coeffs; 
    % Reconstruct using patches' dictionary

    % Add low frequencies to each reconstructed patch
    patches = patches + collect(conf, lores, conf.scale, {});

    % Combine all patches into one image
    img_size = size(lores{1});
    grid = sampling_grid(img_size, conf.window, conf.overlap, conf.border, conf.scale);
    result = overlap_add(patches, img_size, grid);
    % TV desnoise
    if TV_config.TV_enable
        result = TV_2d_denoiser(result, TV_config.TV_lambda, TV_config.TV_maxIter);
    end
    
    % histogram matching
    result = result / mean(result(:)) * mean_lores;
    out_stack(:, :, test_depth) = result;
%     saveastiff(im2uint16(result / max(result(:))), sprintf('depth_%d.tiff', test_depth))
%     fprintf('.');
%     figure(101), imshow(result, [])
end
end

function imgs = modcrop(imgs, modulo)
% imgs as a cell array, since different image might have different size.
for i = 1:numel(imgs)
    sz = size(imgs{i});
    sz = sz - mod(sz, modulo);
    imgs{i} = imgs{i}(1:sz(1), 1:sz(2));
end
end
