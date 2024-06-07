clear all
close all
%/m/nbe/scratch/braindata/bacham1/data_analysis/Analysis//ISC_and_GLM/former_analysis/analysis_230914/dilemma/isc_models.mat

data=load("S:/nbe/braindata/bacham1/data_analysis/Analysis//ISC_and_GLM/former_analysis/analysis_230914/dilemma/isc_models.mat")
%/m/nbe/scratch/braindata/bacham1/data_analysis/Analysis//ISC_and_GLM/former_analysis/analysis_230914/dilemma/isc_input.mat
%/m/nbe/scratch/braindata/bacham1/data_analysis//Articles_niis_and_figs/sisters_study/DATA_FOR_FINAL_VERSION/analysis_2021/svg_figures/Fig_5_eye/Fig_5_eye_analysis_fall2015_basis/full_eISC.mat

mask=triu(ones(120),1);
%mask(find(data.adj_rel==4))=0;
ids=find(mask);
maskperms=zeros(120);
maskperms(ids)=ids;
maskperms=maskperms+maskperms';

for i=1:300
    pe=randperm(120);
    temp=maskperms(pe,pe);
    permIDs(:,i)=temp(ids);
end


%% processing adj_rel
% 1 = acqu
%2 = friend
%3 = sister
%4 = self


temp=data.adj_rel(ids);
model_rel=double([temp==4 temp==3 temp==2 temp==1]); % dummy variables
%model_rel=double([temp==3 temp==2 temp==1]); % dummy variables

model_rel_labels={
'Self'
'Sisters'
'Friends'
'Acquantainces'
};

%% processing adj_pov
% 1 = both have pov of anna, different person
% 2 = both have pov of anna, same person
% 3 = both have pov of kate, different person
% 4  = both have pov of kate, same person

temp=data.adj_pov(ids);
anna=[double(temp==1) + double(temp==2)];
kate=[double(temp==3) + double(temp==4)];
model_pov=[anna kate];
model_pov_labels={
'Anna persp.'
'Kate persp.'
};

%% processing adj_ado
%1 = pairs of non-adopted non same person
%2 = non ado, same
%3 = ado, non same
%4 = ado, same

temp=data.adj_ado(ids);
nonado=[double(temp==1) + double(temp==2)];
ado=[double(temp==3) + double(temp==4)];
model_ado=[nonado ado];
model_ado_labels={
'Non adopt. persp.'
'Adopt. persp.'
};

%% processing eisc
eiscdata=load('S:/nbe/braindata/bacham1/data_analysis//Articles_niis_and_figs/sisters_study/DATA_FOR_FINAL_VERSION/analysis_2021/svg_figures/Fig_5_eye/Fig_5_eye_analysis_fall2015_basis/full_eISC.mat');
temp=eiscdata.full_eISC(ids);
temp(find(temp==0))=NaN;
dummyeisc=double(isnan(temp));
% one approach is to replace nans with medians but take care of this with
% an extra dummy variable
temp(find(isnan(temp)))=nanmedian(temp);
model_eisc=[temp dummyeisc];
model_eisc_labels={
    'eISC'
    'Missing eyetracking data'
};


%% combine models and run regression
model=[model_rel model_pov model_ado model_eisc];

% we do not want to look at the data from the self correlations
poi1=find(model(:,1)==0); % pairs of interest
% we do not want to look at the data with eISC if it has NaNs
poi2=find(~isnan(model(:,end))); % pairs of interest
poi=intersect(poi1,poi2);

model_labels=[
    model_rel_labels
    model_pov_labels
    model_ado_labels
    model_eisc_labels
];


model(:,1) = []; % drop the self
model_labels(1)=[];

allB=[];
allstats=[];
allBsurro=[];
for r=1:264;
    disp(num2str(r))
    load(['isc\iscroi_' num2str(r) '.mat' ]);
    temp=iscmat(ids);
    y=temp(poi);
    [B,BINT,R,RINT,STATS] = regress(zscore(y),[zscore(model(poi,:)) ones(length(poi),1)]);
    %idsp=permIDs(:,r);
    allB(r,:)=B;
    allstats(r,:)=STATS;
    for rrr=1:20
        Bsurro=regress(zscore(y(randperm(length(y)))),[zscore(model(poi,:)) ones(length(poi),1)]);
        allBsurro=[allBsurro;Bsurro'];
    end
end

%% calculate p values from the surrogate null distribution

for reg=1:(size(model,2)-1)
    for r=1:264
        temp=length(find(allB(r,reg)>allBsurro(:,reg)))/size(allBsurro,1); % percentage of strong ones compared to surrogate one
        if(reg==8)
            temp=length(find((allB(r,reg)+allB(r,reg+1))>(allBsurro(:,reg)+allBsurro(:,reg+1))))/size(allBsurro,1); % percentage of strong ones compared to surrogate one
        end
        pvals(r,reg)=1-temp; %right tail p value
    end
    qvals(:,reg)=mafdr(pvals(:,reg),'BHFDR','true');
end

error('stop')

%% net results

load rois_Power264_v2.mat
subnet_ids=[];
for r=1:length(rois)
    subnet_ids(r,1)=rois(r).power_id;
    if(subnet_ids(r,1)>0);
        subnet_labels{subnet_ids(r,1),1}=strrep(rois(r).groupLabel,' ','');
    end
end
subids=[1  3 4 5  7 8 9 10 11 12];
net_labels=subnet_labels(subids);
subids=[0 subids];
net_labels=[
    {'All other ROIs'}
    net_labels];

for subnetblock=1:length(subids)
    thissubid=subids(subnetblock);
    if(thissubid==0)
        roisss=1:264;
    else
        roisss=find(subnet_ids==thissubid);
    end
    for reg=1:8
        allsig(subnetblock,reg)=sum(qvals(roisss,reg)<0.05);
        allsigscaled(subnetblock,reg)=sum(qvals(roisss,reg)<0.05)/length(roisss);
        allmeans(subnetblock,reg)=median(allB(roisss,reg));
        if(reg==8)
            allmeans(subnetblock,reg)=median(allB(roisss,reg)-allB(roisss,reg+1));
        end
    end
end
figure(100)
%subplot(1,2,1)
map=cbrewer('seq','Reds',9);
map=[1 1 1;map];
imagesc(allsigscaled,[0 1])
set(gca,'YTickLabel',net_labels)
set(gca,'XTick',[1:8])
set(gca,'XTickLabel',model_labels(1:8))
set(gca,'XTickLabelRotation',45)

colormap(map)

h=colorbar
ylabel(h,'Percentage of significant regions of interest')

%subplot(1,2,2)
%imagesc(allmeans)
%set(gca,'YTickLabel',net_labels)
%set(gca,'XTick',[1:8])
%set(gca,'XTickLabel',model_labels(1:8))
%set(gca,'XTickLabelRotation',45)%
%
%colormap(map)
%colorbar
%%

close all
Nhist=30
qualmap=cbrewer('qual','Set1',9);
for subnetblock=1:length(subids)
    figure(subnetblock)
    thissubid=subids(subnetblock);
    if(thissubid==0)
        roisss=1:264;
    else
        roisss=find(subnet_ids==thissubid);
    end
    SOFF=10*(subnetblock-1);
    for n=8:-1:1
        
    if(n<8)
        [f x]=ksdensity(allB(roisss,n),-0.2:0.001:0.2);
        allmeans(subnetblock,n)=mean(allB(roisss,n));
    else
        [f x]=ksdensity(allB(roisss,8)+allB(roisss,9),-0.2:0.001:0.2);
        allmeans(subnetblock,n)=mean(allB(roisss,8)+allB(roisss,9));
    end
    hold on
    f=f/10;
    %plot(x,f+n,'LineWidth',1,'Color','white')
    ph=patch(x,f+n-SOFF,1);
    ph.FaceColor=qualmap(n,:)
    ph.FaceAlpha=0.75;
    ph.EdgeColor=[1 1 1]
    %text(-0.1,n+.2-SOFF,model_labels{n})
    end
end

%hf=histogram(allB(:,8),Nhist)
%hf.FaceColor=qualmap(9,:);
set(gca,'YTick',[])
grid on

figure(100)


error('stop')
%% plot
load rois_Power264_v2.mat

subnet_ids=[];
for r=1:length(rois)
    subnet_ids(r,1)=rois(r).power_id;
    if(subnet_ids(r,1)>0);
        subnet_labels{subnet_ids(r,1),1}=strrep(rois(r).groupLabel,' ','');
    end
end
subids=[1  3 4 5  7 8 9 10 11 12];
[aaa bbb]=sort(subnet_ids);
figure(100)
imagesc(allB(bbb,:)')
blocks=find(diff(aaa));
hold on
for b=1:length(blocks)
    %plot([0 4]+.5,[blocks(b) blocks(b)],'k')
    plot([blocks(b) blocks(b)],[0 4]+.5,'k','LineWidth',2)
end

    

