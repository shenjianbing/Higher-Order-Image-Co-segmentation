%----------------------------------------%
function meanFuture = meanFuture(Future,mask)
     meanFuture = [mean(Future(mask,1)) mean(Future(mask,2)) mean(Future(mask,3)) mean(Future(mask,4))];
   