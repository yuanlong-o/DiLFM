clc, clear
close all

addpath('src')
addpath(genpath('ext'))
%% this file is used to demo the DiLFM technique
%  last update: 3/11/2021. YZ

train_data_path = 'data\\train_bead_sample';
test_data_path = 'data\\test_bead_sample';

input_psf_path = 'PSF';
psf_file_name = 'PSF_M_10_NA_0.5_d_750.0_from_-100_to_100_zspac_10_Nnum_7_OSR_3.mat';
%% run RL for training data
gen_RL_capture(input_psf_path, psf_file_name, train_data_path) % training

%% run dictionary low- and high-resolution patch training
% parameters
threshold = 0.01; % ignore blank sample

% dictioanry atom setting
conf.level = 1; % # of scale-ups to perform, what does this me an?
conf.window = [10 10]; % low-res. window size
conf.border = [1 1]; % border of the image (to ignore)
conf.overlap = [9 9]; % partial overlap (for faster training)

% KSVD parameter
ksvd_conf.iternum = 30; % k-svd iteration number
ksvd_conf.memusage = 'normal'; 
ksvd_conf.dictsize = 1500; 
ksvd_conf.Tdata = 5; % maximal sparsity

high_res_dir = sprintf('%s\\high_res_volume_crop', train_data_path);
low_res_dir = sprintf('%s\\LMF_recon', train_data_path);
lfm_dictionary_training(input_psf_path, psf_file_name, threshold, conf, ksvd_conf, ...
                                high_res_dir, low_res_dir)
                            

%% run RL for test data
gen_RL_capture(input_psf_path, psf_file_name, test_data_path) % training
lores_stack = loadtiff(sprintf('%s\\LMF_recon\\sample_1.tif', test_data_path));
lores_stack = double(lores_stack );
%% run scale up
clc
TV_config.TV_enable = false;
TV_config.TV_lambda= 1;
TV_config.TV_maxIter = 10;
dictionary_name = 'data\\Learned_dictionary_size_1500_Tdata_5_peak_4_cut_0.01_overlap_9';
out_stack = lfm_dictionary_test(input_psf_path, psf_file_name, TV_config, dictionary_name, lores_stack, ksvd_conf.Tdata);

if exist('.\\output\\beads.tif' )
    delete ('.\\output\\beads.tif' );
end
saveastiff(im2uint16(out_stack / max(out_stack(:))), '.\\output\\beads.tif');
%% plot
figure, subplot(1, 2, 1), imshow(max(lores_stack, [], 3), []), title('RL')
subplot(1, 2, 2), imshow(max(out_stack, [], 3), []), title('DiLFM')
%% reconstruction
function imgs = modcrop(imgs, modulo)
% imgs as a cell array, since different image might have different size.
for i = 1:numel(imgs)
    sz = size(imgs{i});
    sz = sz - mod(sz, modulo);
    imgs{i} = imgs{i}(1:sz(1), 1:sz(2));
end
end