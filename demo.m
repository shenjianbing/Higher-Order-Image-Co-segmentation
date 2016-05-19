clear all;
close all;

addpath( genpath( '.' ) );

imageFiles = dir(fullfile('./Datasets/image/'));
imageFilesNUM = length(imageFiles)-2;
for imgFilesnum = 1:imageFilesNUM

    %% paths
    img_path = ['Datasets/image/' imageFiles(imgFilesnum+2).name '/image/'];
    scribble_path = ['Datasets/image/' imageFiles(imgFilesnum+2).name '/scribble/'];
    out_path = ['./result/' imageFiles(imgFilesnum+2).name '/'];

    mkdir([out_path 'result/']);
    mkdir([out_path 'regions/']);
    %% parameters
    full_connect = 1;
    lambda = 1;
    hs = 35; hr = 35;  
    M  = 30;

    colors_Bp = [];
    colors_Fp = [];
    lab_colors_Bp = [];
    lab_colors_Fp = [];

    imgFiles = dir(fullfile(img_path ,'*.bmp'));
    LengthimgFiles = length(imgFiles);
    for k = 1:LengthimgFiles;                    
        if( ~exist( fullfile( out_path, 'regions', [imgFiles(k).name(1:end-4), '.mat'] ), 'file' ) )
            img_name = strcat(img_path,imgFiles(k).name);
            img = imread(img_name);
            lab_img = colorspace('Lab<-', img);
            [segs, labels, seg, colors_s, lab_colors_s, edges_s] = msseg(double(img),reshape(lab_img, size(img,1)*size(img,2), size(img,3)),hs,hr,M,full_connect);  
            filename = fullfile( out_path, 'regions', [imgFiles(k).name(1:end-4), '.mat'] );
            save( filename, 'labels', 'edges_s', 'lab_colors_s','-v7.3' );
        end
    end
    
    scribbleFiles = dir(fullfile(scribble_path,'*.bmp'));
    LengthscribbleFiles = length(scribbleFiles);

    %% reading scribbles
    for k = 1:LengthscribbleFiles;
        label_img_name = strcat(scribble_path,scribbleFiles(k).name);%labeling image
        img_name = strcat(img_path,scribbleFiles(k).name);%original image
        img = imread(img_name);
        lab_img = colorspace('Lab<-', img);
        [lines] = seed_generation(label_img_name);
    
        overseg = load( fullfile( out_path, 'regions', [scribbleFiles(k).name(1:end-4), '.mat']));
        labels = overseg.labels;
        lab_colors_s = overseg.lab_colors_s;
    
        labelf = unique(labels(logical(lines(:,1))));%foreground regions
        lab_colors_Fp = [lab_colors_Fp;lab_colors_s(labelf,:)];        
        labelb = unique(labels(logical(lines(:,2))));%background regions 
        lab_colors_Bp = [lab_colors_Bp;lab_colors_s(labelb,:)];      
    end

    %% higher order image cosegmentation
    for k = 1:LengthimgFiles;
        img_name = strcat(img_path,imgFiles(k).name);
        img = imread(img_name);
        fprintf('%s£ºlikelihood estimation.\n',imgFiles(k).name);
        [Posteriors Q Qs n_X n_Y] = likelihood_estimation(lambda,  imgFiles(k).name, img_path, out_path, lab_colors_Fp, lab_colors_Bp);
        imwrite(reshape(Q(:,1),size(img,1),size(img,2)), [out_path 'result/' imgFiles(k).name(1:end-4) '_pri.bmp']);
        fprintf('%s£ºhigher order cosegmentation.\n',imgFiles(k).name);
        higher_order_coseg(imgFiles(k).name,Q,Qs,lab_colors_Fp,lab_colors_Bp,img_path,out_path);  
    end
    
end






