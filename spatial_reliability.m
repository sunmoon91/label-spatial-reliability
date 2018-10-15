function [sr] = spatial_reliability(hard_label_idx, ma, ss)
%calculate the label reliabiltiy
[idx, idy, idz] = ind2sub(size(img), find(ma == 1));
ind = [idx, idy, idz];
[~, ~, ol] = compute_pat(ind,double(hard_label_idx),double(hard_label_idx), ss);
neighbor=compute_patind(ind,double(ma),ss);

sr=zeros(size(ma));
for k=1:size(ind,1)
    a=neighbor(k,:);
    Patch_t=ol(a(a~=0));
    t=ol(k);
    sr(ind(k,1),ind(k,2),ind(k,3))=size((find(Patch_t==t)),1)/size(Patch_t,1);
end
sr(find(ma==0))=-1;
end

