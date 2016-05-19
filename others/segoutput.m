function [imgMasks,segOutline,imgMarkup]=segoutput(img,solution,bi)
%Function [imgMasks,segOutline,imgMarkup]=segoutput(img,solution) 
%generates output for a Cartesian image with dimensions X/Y and a solution 
%conisiting of integer valued nodes indicating membership in a segment
%
%Inputs: img - Original image
%        solution - A 1xN vector assigning an integer to each node
%           indicating its membership in a segment
%
%Outputs: imgMasks - An image where every pixel is assigned an integer 
%           such that pixels sharing numbers belong to the same segment
%         segOutline - A white background with black lines indicating the 
%           segments borders
%         imgMarkup - The same image as the inputs with the red channel 
%           set to 1 along the borders of segments
%
%
%5/23/03 - Leo Grady

% Copyright (C) 2002, 2003 Leo Grady <lgrady@cns.bu.edu>
%   Computer Vision and Computational Neuroscience Lab
%   Department of Cognitive and Neural Systems
%   Boston University
%   Boston, MA  02215
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
% Date - $Id: segoutput.m,v 1.3 2003/08/21 17:29:29 lgrady Exp $
%========================================================================%
if nargin < 3,
    bi = 0;
end;

%Inputs
[X Y Z]=size(img);
nlabels = max(solution(:));
bg_area = find(solution(:)==nlabels);

%Build outputs
imgMasks=reshape(solution,X,Y);

%Outline segments
imgSeg=imgMasks;
if X*Y == 1
    [fx,fy]=deal([]);
elseif X == 1
    fx=[];
    fy=gradient(imgSeg);
elseif Y == 1
    fx=gradient(imgSeg);
    fy=[];    
else
    [fx,fy]=gradient(imgSeg);
end
[xcont_i xcont_j]=find(fx);
[ycont_i ycont_j]=find(fy);
xcont = find(fx);
ycont = find(fy);

segOutline=ones(X,Y);
segOutline(xcont)=0;
segOutline(ycont)=0;

imgMarkup=img(:,:,1); 
if bi == 1, imgMarkup(bg_area)=imgMarkup(bg_area)/3; end;
val1 = 1;
imgMarkup(xcont)=val1;
if xcont_i+1 <= X, imgMarkup(xcont_i+1+X*(xcont_j-1))=val1; end;
if xcont_i-1  > 0, imgMarkup(xcont_i-1+X*(xcont_j-1))=val1; end;
if xcont_j-1  > 0, imgMarkup(xcont_i+X*(xcont_j-2))  =val1; end;
if xcont_j+1 <= Y, imgMarkup(xcont_i+X*(xcont_j))    =val1; end;

if xcont_i+1 <= X, if xcont_j-1  > 0, imgMarkup(xcont_i+1+X*(xcont_j-2))  =val1; end; end;
if xcont_i+1 <= X, if xcont_j+1 <= Y, imgMarkup(xcont_i+1+X*(xcont_j))    =val1; end; end;
if xcont_i-1  > 0, if xcont_j-1  > 0, imgMarkup(xcont_i-1+X*(xcont_j-2))  =val1; end; end;
if xcont_i-1  > 0, if xcont_j+1 <= Y, imgMarkup(xcont_i-1+X*(xcont_j))    =val1; end; end;

imgMarkup(ycont)=val1; 
if ycont_i+1 <= X, imgMarkup(ycont_i+1+X*(ycont_j-1))=val1; end;
if ycont_i-1  > 0, imgMarkup(ycont_i-1+X*(ycont_j-1))=val1; end;
if ycont_j-1  > 0, imgMarkup(ycont_i+X*(ycont_j-2))  =val1; end;
if ycont_j+1 <= Y, imgMarkup(ycont_i+X*(ycont_j))    =val1; end;

if ycont_i+1 <= X, if ycont_j-1  > 0, imgMarkup(ycont_i+1+X*(ycont_j-2))  =val1; end; end;
if ycont_i+1 <= X, if ycont_j+1 <= Y, imgMarkup(ycont_i+1+X*(ycont_j))    =val1; end; end;
if ycont_i-1  > 0, if ycont_j-1  > 0, imgMarkup(ycont_i-1+X*(ycont_j-2))  =val1; end; end;
if ycont_i-1  > 0, if ycont_j+1 <= Y, imgMarkup(ycont_i-1+X*(ycont_j))    =val1; end; end;

if Z == 3
    imgTmp2=img(:,:,2); val2 = 1;
    if bi == 1, imgTmp2(bg_area)=imgTmp2(bg_area)/3; end;
    imgTmp2(xcont)=val2;
    if xcont_i+1 <= X, imgTmp2(xcont_i+1+X*(xcont_j-1))=val2; end;
    if xcont_i-1  > 0, imgTmp2(xcont_i-1+X*(xcont_j-1))=val2; end;
    if xcont_j-1  > 0, imgTmp2(xcont_i+X*(xcont_j-2))  =val2; end;
    if xcont_j+1 <= Y, imgTmp2(xcont_i+X*(xcont_j))    =val2; end;
    
    if xcont_i+1 <= X, if xcont_j-1  > 0, imgTmp2(xcont_i+1+X*(xcont_j-2))  =val2; end; end;
    if xcont_i+1 <= X, if xcont_j+1 <= Y, imgTmp2(xcont_i+1+X*(xcont_j))    =val2; end; end;
    if xcont_i-1  > 0, if xcont_j-1  > 0, imgTmp2(xcont_i-1+X*(xcont_j-2))  =val2; end; end;
    if xcont_i-1  > 0, if xcont_j+1 <= Y, imgTmp2(xcont_i-1+X*(xcont_j))    =val2; end; end;

    imgTmp2(ycont)=val2;
    if ycont_i+1 <= X, imgTmp2(ycont_i+1+X*(ycont_j-1))=val2; end;
    if ycont_i-1  > 0, imgTmp2(ycont_i-1+X*(ycont_j-1))=val2; end;
    if ycont_j-1  > 0, imgTmp2(ycont_i+X*(ycont_j-2))  =val2; end;
    if ycont_j+1 <= Y, imgTmp2(ycont_i+X*(ycont_j))    =val2; end;
    
    if ycont_i+1 <= X, if ycont_j-1  > 0, imgTmp2(ycont_i+1+X*(ycont_j-2))  =val2; end; end;
    if ycont_i+1 <= X, if ycont_j+1 <= Y, imgTmp2(ycont_i+1+X*(ycont_j))    =val2; end; end;
    if ycont_i-1  > 0, if ycont_j-1  > 0, imgTmp2(ycont_i-1+X*(ycont_j-2))  =val2; end; end;
    if ycont_i-1  > 0, if ycont_j+1 <= Y, imgTmp2(ycont_i-1+X*(ycont_j))    =val2; end; end;
    imgMarkup(:,:,2)=imgTmp2;
    
    imgTmp3=img(:,:,3); val3 = 0;
    if bi == 1, imgTmp3(bg_area)=imgTmp3(bg_area)/3; end;
    imgTmp3(xcont)=val3;
    if xcont_i+1 <= X, imgTmp3(xcont_i+1+X*(xcont_j-1))=val3; end;
    if xcont_i-1  > 0, imgTmp3(xcont_i-1+X*(xcont_j-1))=val3; end;
    if xcont_j-1  > 0, imgTmp3(xcont_i+X*(xcont_j-2))  =val3; end;
    if xcont_j+1 <= Y, imgTmp3(xcont_i+X*(xcont_j))    =val3; end;
    
    if xcont_i+1 <= X, if xcont_j-1  > 0, imgTmp3(xcont_i+1+X*(xcont_j-2))  =val3; end; end;
    if xcont_i+1 <= X, if xcont_j+1 <= Y, imgTmp3(xcont_i+1+X*(xcont_j))    =val3; end; end;
    if xcont_i-1  > 0, if xcont_j-1  > 0, imgTmp3(xcont_i-1+X*(xcont_j-2))  =val3; end; end;
    if xcont_i-1  > 0, if xcont_j+1 <= Y, imgTmp3(xcont_i-1+X*(xcont_j))    =val3; end; end;
    
    imgTmp3(ycont)=val3;
    if ycont_i+1 <= X, imgTmp3(ycont_i+1+X*(ycont_j-1))=val3; end;
    if ycont_i-1  > 0, imgTmp3(ycont_i-1+X*(ycont_j-1))=val3; end;
    if ycont_j-1  > 0, imgTmp3(ycont_i+X*(ycont_j-2))  =val3; end;
    if ycont_j+1 <= Y, imgTmp3(ycont_i+X*(ycont_j))    =val3; end;
    
    if ycont_i+1 <= X, if ycont_j-1  > 0, imgTmp3(ycont_i+1+X*(ycont_j-2))  =val3; end; end;
    if ycont_i+1 <= X, if ycont_j+1 <= Y, imgTmp3(ycont_i+1+X*(ycont_j))    =val3; end; end;
    if ycont_i-1  > 0, if ycont_j-1  > 0, imgTmp3(ycont_i-1+X*(ycont_j-2))  =val3; end; end;
    if ycont_i-1  > 0, if ycont_j+1 <= Y, imgTmp3(ycont_i-1+X*(ycont_j))    =val3; end; end;
    imgMarkup(:,:,3)=imgTmp3;
else
    imgTmp1=img(:,:,1);
    imgTmp1(xcont)=0;
    imgTmp1(ycont)=0;
    imgMarkup(:,:,2)=imgTmp1;
    imgMarkup(:,:,3)=imgTmp1;
end
