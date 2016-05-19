
function [W n_Y] = getW_multi(imgValsY,edgesY,beta)

n_Y = size(imgValsY,1);   
imgVals     = imgValsY;
total_edges = edgesY;
n_stack =  n_Y;

weights=makeweights(double(total_edges),imgVals,beta);
W=adjacency(double(total_edges),weights,n_stack);

