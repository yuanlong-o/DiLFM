function [im_h] = backprojection(im_h, im_l, maxIter)
%% for global consistence. See "Image Super-Resolution as Sparse Representation of Raw Image Patches"
%  commented by YZ. last update: 1/26/2020.

[row_l, col_l] = size(im_l);
[row_h, col_h] = size(im_h);

% gaussian blur, aka backprojection filter.
p = fspecial('gaussian', 5, 1);
p = p.^2;
p = p./sum(p(:));


im_l = double(im_l);
im_h = double(im_h);

for ii = 1:maxIter
    % down sampling
    im_l_s = imresize(im_h, [row_l, col_l], 'bicubic');
    
    % calculate difference
    im_diff = im_l - im_l_s;
    
    % upsample again
    im_diff = imresize(im_diff, [row_h, col_h], 'bicubic');
    
    % blur and sum?
    im_h = im_h + conv2(im_diff, p, 'same');
end
    