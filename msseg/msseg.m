% Performing mean_shift image segmentation using EDISON code implementation
% of Comaniciu's paper with a MEX wrapper from Shai Bagon. links at bottom
% of help
%
% Usage:
%   [S L] = msseg(I,hs,hr,M)
%    
% Inputs:
%   I  - original image in RGB or grayscale
%   hs - spatial bandwith for mean shift analysis
%   hr - range bandwidth for mean shift analysis
%   M  - minimum size of final output regions
%
% Outputs:
%   S  - segmented image
%   L  - resulting label map
%
% Links:
% Comaniciu's Paper
%  http://www.caip.rutgers.edu/riul/research/papers/abstract/mnshft.html
% EDISON code
%  http://www.caip.rutgers.edu/riul/research/code/EDISON/index.html
% Shai's mex wrapper code
%  http://www.wisdom.weizmann.ac.il/~bagon/matlab.html
%
% Author:
%  This file and re-wrapping by Shawn Lankton (www.shawnlankton.com)
%  Nov. 2007
%------------------------------------------------------------------------

function [S L seg seg_vals seg_lab_vals seg_edges] = msseg(img,lab_vals,hs,hr,M,full)
  gray = 0;
  if(size(img,3)==1)
    gray = 1;
    I = repmat(img,[1 1 3]);%将矩阵 A 复制 m×n×p 块
  end
  
  if(nargin < 5)
    hs = 10; hr = 7; M = 30;
  end
  
  if(nargin < 6)
    full = 0;
  end
    
  [fimg labels modes regsize grad conf] = edison_wrapper(img,@RGB2Luv,...
      'SpatialBandWidth',hs,'RangeBandWidth',hr,...
      'MinimumRegionArea',M,'speedup',3);
%   fimage  - the result in feature space
%   labels  - labels of regions [if steps==2] 对原图像的过分割标定,对原图像每个点标记
%   modes   - list of all modes [if steps==2]
%   regSize - size, in pixels, of each region [if steps==2]
%   grad    - gradient map      [if steps==2 and synergistic]
%   conf    - confidence map    [if steps==2 and synergistic]
  
  S = fimg; %Luv2RGB(fimg); 
  L = labels + 1; 

  if(gray == 1)
    S = rgb2gray(S);%将真彩色图像转换为灰度图像。
  end
  
  [X,Y,Z] = size(img); nseg = max(L(:)); %有多少个过分割的区域
  vals = reshape(img,X*Y,Z);%重新调整矩阵的行数、列数、维数
  
  if full == 1,
      [x y] = meshgrid(1:nseg,1:nseg);
      %[X,Y] = meshgrid(x,y)
      %将向量x和y定义的区域转换成矩阵X和Y，这两个矩阵可以用来表示mesh和surf的三维空间点以及两个变量的赋值。
      %其中矩阵X的行向量是向量x的简单复制，而矩阵Y的列向量是向量y的简单复制。
      %x=[1 2 3;1 2 3;1 2 3];y=[1 1 1;2 2 2;3 3 3]
      seg_edges = [x(:) y(:)];%区域1和其他区域,区域2和其他区域，。。。。。。
  else
      [points edges]=lattice(X,Y,0);    clear points;
      d_edges = edges(find(L(edges(:,1))~=L(edges(:,2))),:);%找出边缘
      all_seg_edges = [L(d_edges(:,1)) L(d_edges(:,2))]; 
      all_seg_edges = sort(all_seg_edges,2);
  
      tmp = zeros(nseg,nseg);
      tmp(nseg*(all_seg_edges(:,1)-1)+all_seg_edges(:,2)) = 1;
      [edges_x edges_y] = find(tmp==1); 
      seg_edges = [edges_x edges_y];
  end;

  seg_vals = zeros(nseg,Z);%z=3或1
  seg_lab_vals = zeros(nseg,size(lab_vals,2));%size(lab_vals,2)=X*Y,lab_vals图像在lab空间
  for i=1:nseg
    seg{i} = find(L(:)==i);%区域i
    seg_vals(i,:) = mean(vals(seg{i},:));%vals = reshape(img,X*Y,Z);
    seg_lab_vals(i,:) = mean(lab_vals(seg{i},:));%lab空间
  end;
  
  
