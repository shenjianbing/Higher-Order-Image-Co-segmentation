%----------------------------------------%
function dis = disFuture(p,q)
    dis = sqrt(sum((p(:)-q(:)).^2));
   