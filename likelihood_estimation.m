function [Posteriors Qp Qs n_X n_Y ] = likelihood_estimation(lambda, name, img_path, out_path,lab_colors_Fp,lab_colors_Bp)

beta = 60;

img = imread([img_path name]);
n_X = size(img,1)*size(img,2);

overseg = load( fullfile( out_path, 'regions', [name(1:end-4), '.mat']));
labels = overseg.labels;
lab_colors = overseg.lab_colors_s;
edges_s = overseg.edges_s;
%% make affinity matrix 
[W n_Y] = getW_multi(lab_colors,edges_s,beta);
iD = diag(sparse(1./sum(W')));    
PY = iD*W;   
clear iD;

    u = zeros(size(lab_colors,1),2);
    m = zeros(size(lab_colors,1),size(lab_colors_Fp,1));
    for i = 1:size(lab_colors_Fp,1)
        t = lab_colors-repmat(lab_colors_Fp(i,:),size(lab_colors,1),1);
        m(:,i) = sum(t.^2,2); 
    end
    u(:,1)=min(m,[],2)+0.001;

    m = zeros(size(lab_colors,1),size(lab_colors_Bp,1));
    for i = 1:size(lab_colors_Bp,1)
        t = lab_colors-repmat(lab_colors_Bp(i,:),size(lab_colors,1),1);
        m(:,i) = sum(t.^2,2);
    end
    u(:,2)=min(m,[],2)+0.001;
    clear m;
    u_vec = sum(u,2);
    u_vec = u(:,2)./u_vec;

idx(1,1:sum(n_Y))=1;  
I=diag(sparse(idx)); 
p = (I - PY +lambda*I)\(lambda*I*u_vec);
p = [p 1-p];
Q = zeros(n_X,2);
t = zeros(n_X,2);
labels=reshape(labels,n_X,1);

for i = 1:n_Y(1)
    t(labels(:)==i,1) = p(i,1);
    t(labels(:)==i,2) = p(i,2);
end
Q = Q+t;
Qp = Q;
Qs = p(1:n_Y(1),:);
Posteriors = [Q;p];
clear iD;

