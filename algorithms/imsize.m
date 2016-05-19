%----------------------------------------%
function varargout = imsize(img)
iz = size(img);
iz = iz(1:2);
if nargout == 1
    varargout{1} = iz;
else
    varargout{1} = iz(1);
    varargout{2} = iz(2);
end