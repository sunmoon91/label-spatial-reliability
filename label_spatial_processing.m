function [refine_label_pro] = label_spatial_processing(image, label, prob, mask, alpha, opts)
%label-spatial reliability-based label fusion
%image-   intensity image name, type 'image.hdr'&'image.img'
%label-   a set of label
%prob-    a set of label probibality map of segmented using the
%         conventional method, type 'prob_label.hdr'&'prob_label.img'
%mask-    a mask for ROIs segmentation
%alpha-   parameters.
%opts-      Optional inputs (default value: opts=[])


img = analyze75read([image,'.hdr']);

ma = analyze75read([mask,'.hdr']);

prob_labels=zeros([size(img), length(label)]);
for i =1:length(label)
    prob_labels(:,:,:,i)= analyze75read([prob,'_',num2str(label(i)),'.hdr']);
end


[~, hard_label_idx]=max(prob_labels,[],4);

label_reli = label_reliability(prob_labels, label, ma);
spatial_reli = spatial_reliability(hard_label_idx, ma, opts.spatial);

%%%%% calcualting reliabilty with different methods
if opts.rFlag == 1
    reliabilty = label_reli;
end
if opts.rFlag == 2
    reliabilty = spatial_reli;
end
if opts.rFlag == 3
    reliabilty = label_reli .* spatial_reli;
end
if opts.rFlag == 4
    reliabilty = max(label_reli, spatial_reli);
end
if opts.rFlag==5
    reliability = min(label_reli, spatial_reli);
end


intervals = opts.intervals;


[idx, idy, idz] = ind2sub(size(img), find(ma == 1));

ind = [idx, idy, idz];

[P, P1, ol] = compute_pat(ind,double(Img),double(hard_label_idx), opts.patch);

neighbor=compute_patind(ind,double(ma),opts.search);

compute_pat(ind,double(Img),double(I2),p);

iter = floor(1/intervals);
flag=zeros(length(ol));

refine_label_pro=prob_labels;


for i = 1:iter
    range1 = iter-i*intervals;
    if i == iter
        range2 = 0;
    else
        range2 = iter-(i+1)*iteravals;
    end
    
    flag(find(reliability>=range1)) = 1;
    id = find(confidence>=range2 &confidence<range1);
    
    
    for k=1:length(id,1)
        kk=id(k);
        y1=P1(:,kk);
        
        a=neighbor(kk);
        
        Patch1=P1(:,a(a~=0));
        Patch=P(:,a(a~=0));
        Patch_t=ol(a(a~=0));
        Patch_r=reliabilty(a(a~=0));
        
        Patch1=Patch1(1:size(P1,1),:);
        Patch = Patch(1:size(P,1) ,:);
        
        uPatch=logical(flag(a(a~=0)));
        
        Patch1=Patch1(1:size(P1,1),uPatch);
        Patch=Patch(1:size(P,1),uPatch);
        Patch_t=Patch_t(uPatch);
        Patch_r=Patch_r(uPatch);
        
        
        %%%%% normalized the patch
        if opts.sFlag == 1
            uPatch=(abs(mean(Patch1,1)-mean(y1))/mean(y1))<0.01;
            Patch1=Patch1(1:size(P,1),uPatch);
            Patch=Patch(1:size(P,1),uPatch);
            Patch_r=Patch_r(uPatch);
            Patch_t=Patch_t(uPatch)';
        end
        
        
        %%%%% normalized the patch
        if opts.nFlag == 1
            A = Patch;
            y = P(:,kk);
        else
            A = Patch1;
            y = y1;
        end
        
        r=Patch_r;
        t=Patch_t;
        
        if isempty(A)||size(A,2)==1
            continue;
        end
        
        
        A=double(A);
        y=double(y);
        
        %%%%% calcualting voting weight with different methods
        
        if opts.mFlag == 1
            tempw=sum((repmat(y,1,size(A,2))-A).^2);
            m=min(tempw);
            m=sqrt(m);
            w=exp(-1*tempw/(m+0.01)).*r';
        end
        
        if opts.mFlag == 2
            Mij=abs(repmat(y,1,size(A,2))-A);
            M=Mij'*Mij;
            M=M.^opts.beta;
            M=M+0.01*eye(size(M));
            M=M*opts.alpha;
            unit_one=ones(size(A,2),1);
            inv_M=inv(M);
            w=(inv_M*unit_one)/(unit_one'*inv_M*unit_one);
        end
        
        if opts.mFlag == 3
            [w, ~]= nnLeastR(A, y, opts.lambda, opts.lasso);
        end
        
        sum_w =  sum(w);
        refine_prob = zeros(lenght(label));
        for l =1:length(label)
            refine_prob = sum(w(find(t==l)))/sum_w;
            refine_label_pro(ind(k,1),ind(k,2),ind(k,3),l) = (alpha * refine_label_pro(ind(k,1),ind(k,2),ind(k,3),l) + (1-alpha) * refine_prob) * 0.5;
        end
    end
end
end

