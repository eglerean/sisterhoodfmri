close all
clear all
T=712;
R=264;
ids=find(triu(ones(R),1));
surro=[];
for t=10:5:(T-10)
    disp(num2str(t))
    load(['isfc_stats/' num2str(t) '.mat']); %isfc_mat
    surro=[surro;isfc_mat(ids)];
end
[fi xi]=ksdensity(surro,'function','cdf','npoints',5000);

load(['isfc_stats/' num2str(0) '.mat']); %isfc_mat
disp('computing pvalues')
for i=1:length(ids)
    out=isfc_mat(ids(i));
    pval_left=interp1([-1 xi 1],[0 fi 1],out);    % trick to avoid NaNs
    pval=1-pval_left;
    %isfc_mat_pvals(i)=min(pval_left,pval);
    isfc_mat_pvals(i)=pval;
end

isfc_mat_qvals=mafdr(isfc_mat_pvals,'BHFDR','true');

save('isfc_results.mat','isfc_mat','isfc_mat_qvals');