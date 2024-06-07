clear all
close all
T=712;
R=264;

%% for the 1st view of the movie, let's compute surrogate ISFC adjacency matrices that we will use for computing the null distribution
for view=1:1
    parfor t=1:T
        isfc_compute(view,t);
    end
end

%% helper function to compute surrogate ISFC adj matrices

function isfc_compute(view,t)
    disp(['Computing view ' num2str(view) ' t ' num2str(t)])
    t=t-1;
    isfc_mat=zeros(264);
    % first we compute a fake ISFC adj matrix for each subject
    for s=1:30
        load(['isfc_data/view' num2str(view) '/this_' num2str(s) '.mat']); % thisTS
        load(['isfc_data/view' num2str(view) '/others_' num2str(s) '.mat']) % othersTS
        % we circular shift the N-1 average subjects used for this subject
        % ISFC
        othersTS=circshift(othersTS,t);
        adj=corr([thisTS othersTS]);
        % the ISFC matrix is the 264x264 block on top right
        subb=atanh(adj(1:264,265:end));
        % we need to make the ISFC matrix symmetric
        this_isfc=(subb+subb')/2; 
        % we add it to the surrogate matrix across all subjects
        isfc_mat=isfc_mat+this_isfc;
    end
    isfc_mat=tanh(isfc_mat/30); % we average across subjects and z-transform invert
    save(['isfc_stats/' num2str(t) '.mat'],'isfc_mat')
end

