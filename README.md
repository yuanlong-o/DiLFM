# DiLFM

Implementation for dictionary lgiht field microscope (DiLFM)

## System requirements
* Software:
  * Matlab-compatible version of Windows or Linux (see https://www.mathworks.com/support/requirements/matlab-system-requirements.html)
  * Matlab R2019b or newer
  * Toolboxes: signal processing, image processing
  * Matlab-compatible CUDA driver (for LFM deconvolution)

## Usage
* To run the DiLFM demo
  1. Generate /src/gen_sphere_sample.m for train and test sample generation
  2. Download example psf from Google drive https://drive.google.com/file/d/1XtifZ0tZmhtqDTLShV8mhZz-lIJXyH9y/view?usp=sharing
  3. Run demo_dictionary_LFM.m and check out_stack varaible for final output

## Instructions
* To train your own DiLFM model
  1. Collect trainning data that is similar with your target data
  2. Run virtual LFM propagation and RL deconvolution to get low-resolution results (this can be done through /src/gen_RL_capture.m)
  3. Run lfm_dictionary_training.m for high- and low-resolution pair training
