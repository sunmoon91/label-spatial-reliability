function [ lr ] = label_reliability(label_prob_map, label, ma)
%calculate the label reliabiltiy
%the label reliabiltiy will be normalized
lr=zeros(size(ma));
for i = 1:lenght(label)
    prob = label_prob_map(:,:,:,i);
    lr = lr + prob.*log2(prob+0.0001);
end

max_lr=max(max(max(lr)));
min_lr=min(min(min(lr)));
lr = (lr-min_lr)./(max_lr-min_lr);
lr(find(ma==0))=-1;
end

