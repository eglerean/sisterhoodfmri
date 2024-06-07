close all
clear all

% remote folder
% addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila//bramila'));
% local folder
addpath(genpath('/nbe/braindata/shared/toolboxes/bramila//bramila'));

map=cbrewer('seq','Reds',9);
map=[1 1 1;map];
load('isfc_results.mat');
load('rois_Power264_v2.mat')
R=264;
ids=find(triu(ones(R),1));

isfc_mat_vec=isfc_mat(ids);
th=min(isfc_mat_vec(find(isfc_mat_qvals<0.05)));
isfc_mat_th=isfc_mat.*(isfc_mat>th);


%imagesc(isfc_mat_th)
%colormap((map))
%colorbar


%% plot the ISFC matrix, rearranged using Power et al subnetworks

% get subnetowkr IDs
subnet_ids=[];
for r=1:length(rois)
    subnet_ids(r,1)=rois(r).power_id;
    if(subnet_ids(r,1)>0);
        subnet_labels{subnet_ids(r,1),1}=strrep(rois(r).groupLabel,' ','');
    end
end
subids=[1  3 4 5  7 8 9 10 11 12];
[aaa bbb]=sort(subnet_ids);
outnet=zeros(R);
outnet=isfc_mat_th;
%outnet=outnet+outnet';
outnet=outnet(bbb,bbb);

% plot
imagesc(outnet,[0 0.5])
blocks=find(diff(aaa));
hold on
for b=1:length(blocks)
    plot([0 R],[blocks(b) blocks(b)],'k')
    plot([blocks(b) blocks(b)],[0 R],'k')
end

axis square
colormap(map)
h=colorbar
ylabel(h,'ISFC')
intervals=[[.5; blocks-.5  ] [blocks+.5; R+.5]]
outlabels=subnet_labels(subids);
outlabels=[{'n/a'}; outlabels];

set(gca,'XTick',mean(intervals,2));
set(gca,'YTick',mean(intervals,2));
set(gca,'XTickLabel',outlabels)
set(gca,'XTickLabelRotation',45)
set(gca,'YTickLabel',outlabels)
%% do a summary based on percentage of links per subnetwork

for ri=1:length(subids)
    for ci=1:length(subids)
        rows=bbb(find(aaa==subids(ri)));
        cols=bbb(find(aaa==subids(ci)));
        subsub=isfc_mat_th(rows,cols);
        TL=length(subsub(:));
        SL=sum(subsub(:)>0);
        sum_stats(ri,ci)=SL/TL;
       
    end
end
figure
h=imagesc(sum_stats,[0 1])
h.AlphaData=1-triu(ones(10),1)
colormap(map)
axis square
h=colorbar
ylabel(h,'Percentage of significant links')
set(gca,'XTick',1:length(sum_stats));
set(gca,'YTick',1:length(sum_stats));

set(gca,'XTickLabel',subnet_labels(subids))
set(gca,'XTickLabelRotation',45)
set(gca,'YTickLabel',subnet_labels(subids))
%% output some insteresting nodes based on node density

ND=sum(isfc_mat);
load rois_Power264_v2.mat
[aa bb]=sort(ND,'descend');
for i=1:50
    disp([num2str(aa(i)) ' ' rois(bb(i)).groupLabel ' ' num2str(rois(bb(i)).power_id)]);
end

power_ids=[];
for r=1:264
    power_ids(r)=rois(r).power_id;
end

%figure
%imagesc(power_ids)





