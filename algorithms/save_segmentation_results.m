function save_segmentation_results(Posteriors, n_X, name, img_path, out_path, bg_sub)

img = im2double(imread([img_path name]));    h = size(img,1);    w = size(img,2);
nlabels = size(Posteriors,2);

if nargin < 6, bg_sub = 1; end;

%% Display the result of the segmentation for pixels
data_path = [out_path 'segments/']; 
mkdir(data_path);
[vals,inds] = max(Posteriors(1:n_X,:)');
mask = reshape(inds,h,w);

    bound_color = [255 255 0];
    % Highlight boundaries
    I = imread([img_path name]);
    P = bwperim(mask > 1);
    I = apply_color(I, P, bound_color);
    
    function I = apply_color(I, mask, color)
        for i=1:3
            select = zeros([h w 3]);
            select(:, :, i) = mask;         
            I(select ~= 0) = color(i);
        end
    end
imwrite(I, [out_path 'segments/' name '_ours1.bmp']);

[imgMasks,segOutline,imgMarkup]=segoutput(img,mask,bg_sub);
%figure; imshow(imgMarkup);  

clear imgMasks segOutline imgMarkup;

disp_colors(1,1:3) = [1,1,1];
disp_colors(2,1:3) = [0,0,0];
disp_colors(3,1:3) = [0,0,1];
disp_colors(4,1:3) = [1,1,0];
disp_colors(5,1:3) = [1,0,1];
disp_colors(6,1:3) = [0,1,1];
disp_colors(7,1:3) = [1,0,0];
disp_colors(8,1:3) = [0,1,0];
data_path = [out_path 'labels/']; mkdir(data_path);
if nlabels <= 8,
    label_img = zeros(h,w,3);
    for nc=1:3,
        tmp = label_img(:,:,nc);
        for i=1:nlabels, tmp(find(mask==i)) = disp_colors(i,nc); end;
        label_img(:,:,nc) = tmp; clear tmp;
    end;
    imwrite(label_img, [out_path 'labels/' name(1:end-4) '.bmp']);
end;
end
