function lfm_dictionary_training(input_psf_path, file_name, threshold, conf, ksvd_conf, ...
                                high_res_dir, low_res_dir)
%  last update: 4/10/2020. YZ

% input_psf_path = '..\\..\\..\\data_from_xb\\psf\\40x water';
% file_name = 'PSF_M_47_NA_1.0_fobj_3.55mm_d_2100.0_from_-30_to_30_zspac_2_Nnum_13_OSR_3.mat';
load(sprintf('%s/%s', input_psf_path, file_name));

%% parameters
% threshold = 0.01;
% 
% 
% conf.level = 1; % # of scale-ups to perform, what does this me an?
% conf.window = [5 5]; % low-res. window size
% conf.border = [1 1]; % border of the image (to ignore)
% conf.overlap = [4 4]; % partial overlap (for faster training)
% 
% % KSVD parameter
% ksvd_conf.iternum = 30; % TBD
% ksvd_conf.memusage = 'normal'; % higher usage doesn't fit...
% ksvd_conf.dictsize = 1500; % TBD
% ksvd_conf.Tdata = 5; % maximal sparsity: TBD


%% define scale, which determine the patch size
peak = round(Nnum / 2);
depth_N = size(H ,5);
scale(1 : ceil(depth_N / 2) - 1) = round(linspace(round(peak  / 2), peak - 2, ceil(depth_N / 2) - 1));
scale(ceil(depth_N / 2) + 1 : depth_N) = round(linspace( peak  - 2, round(peak  / 2), ceil(depth_N / 2) - 1));
scale(ceil(depth_N / 2)) = peak;


%% generate low and high resolution pair
% high_res_dir = 'high_res_volume';
% low_res_dir = 'LMF_recon'; % the file name should be the same

file_pattern = fullfile(high_res_dir, '*.tif');
tiff_file = dir(file_pattern);

high_res_img_cell = cell(size(H, 5), length(tiff_file));
low_res_img_cell = cell(size(H, 5), length(tiff_file));

for k = 1 : length(tiff_file)
    base_filename = tiff_file(k).name;
    high_res_stack = double(loadtiff(sprintf('%s\\%s',high_res_dir, base_filename)))/ 65536;
    low_res_stack = double(loadtiff(sprintf('%s\\%s',low_res_dir, base_filename))) / 65536;
    
    for kk = 1 : size(H ,5)
        high_res_img_cell{kk, k} = high_res_stack(:, :, kk);
        low_res_img_cell{kk, k} = low_res_stack(:, :, kk);
    end

end

output_folder = sprintf('data\\Learned_dictionary_size_%d_Tdata_%d_peak_%d_cut_%g_overlap_%d', ...
    ksvd_conf.dictsize, ksvd_conf.Tdata,peak, threshold, conf.overlap(1));
mkdir(output_folder)
%% learn the dictionary depth by depth
for ii = 1 : depth_N
    ii
    % low res image
    lores =low_res_img_cell(ii, :);
    
	% high res image
    hires = high_res_img_cell(ii, :); % here since the low res image has a crop
  
    % convolutional filters
    conf.scale = scale(ii);
    % High-pass filters for feature extraction (defined for upsampled low-res.)
    O = zeros(1, conf.scale-1);
    G = [1 O -1]; % Gradient
    L = [1 O -2 O 1]/2; % Laplacian
    
    conf.filters = {G, G.', L, L.'}; % 2D versions
    features = collect(conf, lores, conf.scale, conf.filters);
    
    % judging patches and abandon some, this is different from the natural
    patches_buf = collect(conf, hires, conf.scale, {});
    patches_buf_max = max(patches_buf, [], 1);
    valid_ind = patches_buf_max > threshold;
%     sum(valid_ind)
    % patch
    patches = cell(size(hires));
    for i = 1:numel(patches) % Remove low frequencies
        patches{i} = hires{i} - lores{i};
    end

    patches = collect(conf, patches, conf.scale, {}); % no convolutional kernals, so just features
    patches = patches(:, valid_ind);
    features = features(:, valid_ind);
    
    % PCA dimensionality reduction
    C = double(features * features');
    [V, D] = eig(C);
    D = diag(D); % perform PCA on features matrix 
    D = cumsum(D) / sum(D); % normalize is necessary
    k = find(D >= 1e-3, 1); % ignore 0.1% energy
    conf.V_pca = V(:, k:end); % choose the largest eigenvectors' projection
    conf.ksvd_conf = ksvd_conf;
    features_pca = conf.V_pca' * features;

    % Combine into one large training set
%     clear C D V
    ksvd_conf.data = double(features_pca);
%     clear features_pca
    % Training process (will take a while)
    tic;
    fprintf('Training [%d x %d] dictionary on %d vectors using K-SVD\n', ...
        size(ksvd_conf.data, 1), ksvd_conf.dictsize, size(ksvd_conf.data, 2))
    % run KSVD
    [conf.dict_lores, gamma] = ksvd(ksvd_conf); % , ''); 
    toc;
    % X_lores = dict_lores * gamma
    % X_hires = dict_hires * gamma {hopefully}

    fprintf('Computing high-res. dictionary from low-res. dictionary\n');
    % dict_hires = patches / full(gamma); % Takes too much memory...
    patches = double(patches); % Since it is saved in single-precision.
    dict_hires = (patches * gamma') * inv(full(gamma * gamma'));
    conf.dict_hires = double(dict_hires); 
    
    % save
    save(sprintf('%s\\depth_%d.mat', output_folder, ii), 'conf');
end
end
%% utility function
function imgs = modcrop(imgs, modulo)
% imgs as a cell array, since different image might have different size.
for i = 1:numel(imgs)
    sz = size(imgs{i});
    sz = sz - mod(sz, modulo);
    imgs{i} = imgs{i}(1:sz(1), 1:sz(2));
end
end


function imgs = load_images(paths)

imgs = cell(size(paths));
for i = 1:numel(paths)
    X = imread(paths{i});
    if size(X, 3) == 3 % we extract our features from Y channel
        X = rgb2ycbcr(X);
        X = X(:, :, 1);
    end
    X = im2single(X); % to reduce memory usage
    X = X - min(X(:));
    X = X / max(X(:));
    imgs{i} = single(X);
end
end