clear all
close all
T=712;
R=264;

for view=1:1
    parfor t=1:T
        isfc_compute(view,t);
    end
end

function isfc_compute(view,t)
    disp(['Computing view ' num2str(view) ' t ' num2str(t)])
    t=t-1;
    isfc_mat=zeros(264);
    for s=1:30
        load(['isfc_data/view' num2str(view) '/this_' num2str(s) '.mat']); % thisTS
        load(['isfc_data/view' num2str(view) '/others_' num2str(s) '.mat']) % othersTS
        othersTS=circshift(othersTS,t);
        adj=corr([thisTS othersTS]);
        subb=atanh(adj(1:264,265:end));
        this_isfc=(subb+subb')/2; 
        isfc_mat=isfc_mat+this_isfc;
    end
    isfc_mat=tanh(isfc_mat/30);
    save(['isfc_stats/' num2str(t) '.mat'],'isfc_mat')
end

