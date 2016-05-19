%----------------------------------------%
function Future = getFuture(img)
    [nX nY nZ] = size(img);
    Ix2 = zeros(nX,nY);
    Ixy = zeros(nX,nY);
    Iy2 = zeros(nX,nY);  
    I2 = zeros(nX,nY);
    for i = 1:nZ
    g = double(img(:,:,i));
    [Ix,Iy] = gradient(g);
    Ix(Ix(:)==0)=0.01;
    Iy(Iy(:)==0)=0.01;
    Ix2=Ix.*Ix+Ix2;%Ix?
    Ixy=Ix.*Iy+Ixy;%Iy?
    Iy2=Iy.*Iy+Iy2;%Ixy
    g = g/10;
    I2 =g.*g +I2;
    end
    Future = [I2(:) Ix2(:) Iy2(:) Ixy(:)];
   