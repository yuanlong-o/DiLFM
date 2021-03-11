clc, clear
close all

addpath(genpath('..\\ext'))
%% this file is used to generate a 3D sphere sample 
%  last update: 4/26/2020. 

file_path = '..\\psf';
file_name = 'PSF_M_10_NA_0.5_d_750.0_from_-100_to_100_zspac_10_Nnum_7_OSR_3.mat';
load(sprintf('%s/%s', file_path, file_name));


output_dir_train = '..\\data\\train_bead_sample';
train_sample_number = 20;
mkdir(output_dir_train)

output_dir_test = '..\\data\\test_bead_sample';
test_sample_number  = 1;
mkdir(output_dir_test)

%% parameters
size_x = 280; 
size_y = 280;
size_z = length(x3objspace);
if mod(size_x, Nnum)
    error('size_x must be times of Nnum')
end

pixel_sample = pixelPitch / M;

sphere_size = 20e-6;
N_sphere = 20;
z_ratio = 2;

%% generate sample
sphere_in_pixel_lateral = ceil(sphere_size / pixel_sample);
if ~mod(sphere_in_pixel_lateral, 2)
    sphere_in_pixel_lateral = sphere_in_pixel_lateral + 1;
end
sphere_in_pixel_axial = ceil(sphere_size / zspacing * z_ratio);
if ~mod(sphere_in_pixel_axial, 2)
    sphere_in_pixel_axial = sphere_in_pixel_axial + 1;
end

half_size_lateral = (sphere_in_pixel_lateral - 1)/2;
half_size_axial = (sphere_in_pixel_axial - 1)/2; 
neuron_temp = zeros(sphere_in_pixel_lateral, sphere_in_pixel_lateral, sphere_in_pixel_axial);
x_temp = (1 : sphere_in_pixel_lateral) - half_size_lateral - 1;
z_temp = (1 : sphere_in_pixel_axial) - half_size_axial - 1;
[X_temp, Y_temp, Z_temp] = meshgrid(x_temp, x_temp, z_temp);
neuron_temp((X_temp/sphere_in_pixel_lateral * 2).^2+(Y_temp/sphere_in_pixel_lateral* 2).^2 + (Z_temp / sphere_in_pixel_axial* 2).^2 < 1) = 1;

%% train data
for kk = 1 : train_sample_number
    kk
range_xy = 0.8;
range_z = 1.05;
neuron_pos = rand(N_sphere, 3);
neuron_pos(:, 1 : 2) = round(((neuron_pos(:, 1 : 2) - 0.5) * range_xy + 0.5) * size_x);
neuron_pos(:, 3) = round(((neuron_pos(:, 3) - 0.5) * range_z + 0.5) * size_z);
sample = zeros(size_x, size_y, size_z);
for i = 1 : N_sphere
    ind_x = neuron_pos(i, 1);
    ind_y = neuron_pos(i, 2);
    ind_z = neuron_pos(i, 3);
    % consider the volume boundary
    start_x = half_size_lateral;
    end_x = half_size_lateral;
    if ind_x - half_size_lateral < 1
        start_x = ind_x - 1;
    elseif ind_x + half_size_lateral > size_x
        end_x = size_x - ind_x;
    end
    start_y = half_size_lateral;
    end_y = half_size_lateral;
	if ind_y - half_size_lateral < 1
        start_y = ind_y - 1;
    elseif ind_y + half_size_lateral > size_y
        end_y = size_y - ind_y;
    end
    start_z = half_size_axial;
    end_z = half_size_axial;    
	if ind_z - half_size_axial < 1
        start_z = ind_z - 1;
    elseif ind_z + half_size_axial > size_z
        end_z = size_z - ind_z;
    end
    
    sample(ind_x - start_x : ind_x + end_x , ...
        ind_y - start_y: ind_y + end_y, ...
        ind_z - start_z : ind_z + end_z) = ...
    sample(ind_x - start_x : ind_x + end_x, ...
        ind_y - start_y: ind_y + end_y, ...
        ind_z - start_z : ind_z + end_z)+ ...
    neuron_temp(half_size_lateral + 1 - start_x : half_size_lateral + 1  + end_x, ...
                half_size_lateral + 1 - start_y : half_size_lateral + 1  + end_y, ...
                half_size_axial + 1 - start_z : half_size_axial + 1  + end_z);
end
sample = double(sample > 0.1);
% save
saveastiff(im2uint16(sample / max(sample(:))), sprintf('%s\\sample_%d.tif', ...
    output_dir_train, kk))

end


%% test data
for kk = 1 : test_sample_number
    kk
range_xy = 0.8;
range_z = 1.05;
neuron_pos = rand(N_sphere, 3);
neuron_pos(:, 1 : 2) = round(((neuron_pos(:, 1 : 2) - 0.5) * range_xy + 0.5) * size_x);
neuron_pos(:, 3) = round(((neuron_pos(:, 3) - 0.5) * range_z + 0.5) * size_z);
sample = zeros(size_x, size_y, size_z);
for i = 1 : N_sphere
    ind_x = neuron_pos(i, 1);
    ind_y = neuron_pos(i, 2);
    ind_z = neuron_pos(i, 3);
    % consider the volume boundary
    start_x = half_size_lateral;
    end_x = half_size_lateral;
    if ind_x - half_size_lateral < 1
        start_x = ind_x - 1;
    elseif ind_x + half_size_lateral > size_x
        end_x = size_x - ind_x;
    end
    start_y = half_size_lateral;
    end_y = half_size_lateral;
	if ind_y - half_size_lateral < 1
        start_y = ind_y - 1;
    elseif ind_y + half_size_lateral > size_y
        end_y = size_y - ind_y;
    end
    start_z = half_size_axial;
    end_z = half_size_axial;    
	if ind_z - half_size_axial < 1
        start_z = ind_z - 1;
    elseif ind_z + half_size_axial > size_z
        end_z = size_z - ind_z;
    end
    
    sample(ind_x - start_x : ind_x + end_x , ...
        ind_y - start_y: ind_y + end_y, ...
        ind_z - start_z : ind_z + end_z) = ...
    sample(ind_x - start_x : ind_x + end_x, ...
        ind_y - start_y: ind_y + end_y, ...
        ind_z - start_z : ind_z + end_z)+ ...
    neuron_temp(half_size_lateral + 1 - start_x : half_size_lateral + 1  + end_x, ...
                half_size_lateral + 1 - start_y : half_size_lateral + 1  + end_y, ...
                half_size_axial + 1 - start_z : half_size_axial + 1  + end_z);
end
sample = double(sample > 0.1);
% save
saveastiff(im2uint16(sample / max(sample(:))), sprintf('%s\\sample_%d.tif', ...
    output_dir_test, kk))
end
%% with noise


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