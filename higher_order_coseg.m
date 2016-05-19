function higher_order_coseg(name,Posteriors,Qs,lab_colors_Fp,lab_colors_Bp,img_path,out_path ) 
  N = 8; % neighborhood
  lambda = 20;
  alpha = 50;
  label_cnt = 2;%label numer
  
  img_name = [img_path name];

  s.img=imread(img_name);
  [nX,nY,nZ]=size(s.img);

  s.Future = getFuture(s.img);  
  s.weight = zeros(nX*nY,1);

  lab_img = colorspace('Lab<-', s.img);
overseg = load( fullfile( out_path, 'regions', [name(1:end-4), '.mat']));
labels = overseg.labels;
edges_s = overseg.edges_s;
  
  
  pb = Posteriors(1:size(s.img,1)*size(s.img,2),2);
  pf = Posteriors(1:size(s.img,1)*size(s.img,2),1);
  mask = logical(pb./pf>=0.7&pb./pf<=1.3);
    
    lab_img = reshape(lab_img,size(lab_img,1)*size(lab_img,2),size(lab_img,3));
    lab_colors=lab_img(mask,:);
    u = zeros(size(lab_colors,1),2);
    for i = 1:size(lab_colors_Fp,1)
        t = lab_colors-repmat(lab_colors_Fp(i,:),size(lab_colors,1),1);
        m(:,i) = sum(t.^2,2); 
    end
    u(:,1)=min(m,[],2);
    clear m;
    for i = 1:size(lab_colors_Bp,1)
        t = lab_colors-repmat(lab_colors_Bp(i,:),size(lab_colors,1),1);
        m(:,i) = sum(t.^2,2);
    end
    u(:,2)=min(m,[],2);
    clear m;
    u_vec = sum(u,2);
    u_vec = [u(:,2)./u_vec u(:,1)./u_vec];
    
    pb(mask(:))=u_vec(:,2);
    pf(mask(:))=u_vec(:,1);
    l = zeros(nX,nY);
    l(pf(:)>0.5) = 1;
    BTW=bwconncomp(reshape(l,nX,nY),8);
    for i=1:BTW.NumObjects
     if numel(BTW.PixelIdxList{i}) < nX*nY*0.01
         pf(BTW.PixelIdxList{i}(:))=0.3;
         pb(BTW.PixelIdxList{i}(:))=0.7;
    end
    end
    l = zeros(nX,nY);
    l(pf(:)>0.5) = 1;
    BTW=bwconncomp(~l,8);
    for i=1:BTW.NumObjects
        if numel(BTW.PixelIdxList{i}) <  nX*nY*0.01
            pf(BTW.PixelIdxList{i}(:))=0.7;
            pb(BTW.PixelIdxList{i}(:))=0.3;
        end
    end
 
  V = 1 - eye(label_cnt); % uniform label bias (generally ignore this)
  hop = H();
  map = robustpn_mex([nX nY 1], N, label_cnt,@R, @B,hop,V);
  s.map = reshape(map,size(s.img,1),size(s.img,2));
  s.map = bwmorph(s.map,'clean');
  BTW=bwconncomp(s.map,8);
  for i=1:BTW.NumObjects
     if numel(BTW.PixelIdxList{i}) <  nX*nY*0.01
         s.map(BTW.PixelIdxList{i}(:))=0;
     end
  end
  BTW=bwconncomp(~s.map,8);
  for i=1:BTW.NumObjects
     if numel(BTW.PixelIdxList{i}) <  nX*nY*0.01
         s.map(BTW.PixelIdxList{i}(:))=1;
     end
  end

  rimg =reshape(s.img,nX*nY,1,nZ);
  rimg(s.map==0,1)=255;
  rimg(s.map==0,2)=255;
  rimg(s.map==0,3)=255;
  rimg =reshape(rimg,nX,nY,nZ);
  imwrite(rimg, [out_path 'result/' name '_seg.bmp']);
  p  = zeros(size(map,2),2);
  p(s.map(:)==1,1)=1;
  p(s.map(:)==0,2)=1;
  save_segmentation_results(p, nX*nY, name, img_path, [out_path 'result/']);
  clear all;

  function w = R() % regional term
    wb = -log(pb + eps);
    wf = -log(pf + eps);
    w = alpha*[wb';wf'];
  end

  function w = B(e, p, q) % boundary term
    t = zeros(1,size(p,2));
    w = t;
    if e(3) == 16,
       e(3) = 0;
       for i = 1:size(p,2);
         t(i)=disFuture(s.Future(p(i),:),s.Future(q(i),:));
          if s.weight(p(i)) == 0
              s.weight(p(i)) = 0.1;
          end
          w(i) = lambda * exp( -t(i)/100)/ norm(e);
       end;
    else
      for i = 1:size(p,2);
        t(i)=disFuture(s.Future(p(i),:),s.Future(q(i),:));
        s.weight(p(i))=s.weight(p(i))+t(i);
        s.weight(q(i))=s.weight(q(i))+t(i);
      end;
    end;
  end

function hop = H() % higher order term
    mss = labels;
    lrs=imresize(mss,2*imsize(mss),'bilinear');
    nrs=imresize(mss,2*imsize(mss),'nearest');
    d = bwdist(lrs~=nrs);
    w = imresize(.5*d, imsize(mss), 'bilinear');
    w = w+.1;
    st = regionprops(mss,'PixelIdxList');
    gamma(label_cnt+1) = 10; 
    [hop(1:numel(st))] = deal(struct('ind',[],'w',[],...
    'gamma',single(gamma),'Q',.1));
    for hi=1:numel(st)
      hop(hi).ind = st(hi).PixelIdxList;      
      segMean = meanFuture(s.Future,hop(hi).ind);
      dis = zeros(numel(hop(hi).ind),1);
      for i = 1:numel(hop(hi).ind)     
        dis(i)=disFuture(s.Future(hop(hi).ind(i),:),segMean)+1e-3;
      end
      C = numel(st(hi).PixelIdxList);
      dis = dis .* C ./ sum(dis);
      hop(hi).w = w(hop(hi).ind)./dis;
      hop(hi).w = min(hop(hi).w, 2);
      P = sum(hop(hi).w);
      hop(hi).w = single(hop(hi).w .* C ./ P);%Q = sum(m) = C;w = w/ava(w);
      hop(hi).Q = single(.1 * C);
      hop(hi).gamma(end) = single(C);
      RegionPb = sum(pb(hop(hi).ind));
      RegionPf = sum(pf(hop(hi).ind));
      x = min(C,100)/120;
      region = [edges_s(edges_s(:,1)==hi,2);edges_s(edges_s(:,2)==hi,1)];
      rb = C*(x*RegionPb/C+(1-x)*sum(Qs(region,2))/size(region,1));
      rf = C*(x*RegionPf/C+(1-x)*sum(Qs(region,1))/size(region,1));
      hop(hi).gamma(1) = hop(hi).gamma(end)-rb;
      hop(hi).gamma(2) = hop(hi).gamma(end)-rf;
    end
end
end
