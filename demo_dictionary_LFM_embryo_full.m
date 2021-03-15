clc, clear
close all

addpath('src')
addpath(genpath('ext'))
%% this file is used to demo the DiLFM technique with embryo dataset
%  last update: 3/11/2021. YZ

train_data_path = 'data\\train_bio_sample';
test_data_path = 'data\\test_bio_sample';

input_psf_path = 'PSF';
psf_file_name = 'PSF_M_47_NA_1.0_fobj_3.55mm_d_2100.0_from_-30_to_30_zspac_2_Nnum_13_OSR_3.mat';
%% run RL for training data
gen_RL_capture(input_psf_path, psf_file_name, train_data_path) % training
depth_range = [14, 19];
%% run dictionary low- and high-resolution patch training
% parameters
threshold = 0.01; % ignore blank sample

% dictioanry atom setting
conf.level = 1; % # of scale-ups to perform, what does this me an?
conf.window = [5 5]; % low-res. window size
conf.border = [1 1]; % border of the image (to ignore)
conf.overlap = [4 4]; % partial overlap (for faster training)

% KSVD parameter
ksvd_conf.iternum = 30; % k-svd iteration number
ksvd_conf.memusage = 'normal'; 
ksvd_conf.dictsize = 1500; 
ksvd_conf.Tdata = 5; % maximal sparsity

high_res_dir = sprintf('%s\\high_res_volume_crop', train_data_path);
low_res_dir = sprintf('%s\\LMF_recon', train_data_path);
lfm_dictionary_training(input_psf_path, psf_file_name, threshold, conf, ksvd_conf, ...
                                high_res_dir, low_res_dir, depth_range)
                            

%% load test data
lores_stack = loadtiff(sprintf('%s\\recon_40x_embryo_stack.tif', test_data_path));
lores_stack = double(lores_stack );
%% run dictionary
clc
TV_config.TV_enable = false;
TV_config.TV_lambda= 1;
TV_config.TV_maxIter = 10;
dictionary_name = 'data\\Learned_dictionary_size_1500_Tdata_5_peak_7_cut_0.01_overlap_4';
out_stack = lfm_dictionary_test(input_psf_path, psf_file_name, TV_config, dictionary_name, lores_stack, ksvd_conf.Tdata, depth_range);

if exist('.\\output\\embryo.tif' )
    delete ('.\\output\\embryo.tif' );
end
saveastiff(im2uint16(out_stack / max(out_stack(:))), '.\\output\\embryo.tif');
