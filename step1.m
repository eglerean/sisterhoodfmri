clear all
close all 
load /m/nbe/scratch/braindata/bacham1/data_analysis/Analysis/ISC_and_GLM/ISC271114/Test.mat

for n=1:120
    temp=Params.PublicParams.subjectSource{n};
    subjpath{n}=strrep(temp,'data_analysis','data_analysis/Data');
    if isfile(subjpath{n})
        disp("ok!")
    else
        error(subjpath{n})
    end
end

disp('Adding toolbox')
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila//bramila'));

Nsubj=length(subjpath);
load rois_Power264;
mkdir('rois')

%disp('parfor starts')
%parfor s=1:Nsubj
%    roi_extract_helper(s,subjpath)
%end
if(0)
mkdir('isfc_data')
for view=1:4
    
    temp=1:Nsubj;
    view_vec=temp(view:4:end);
    mkdir(['isfc_data/view' num2str(view)])
    for s=1:Nsubj/4
        
        this=view_vec(s);
        disp(['view ' num2str(view) ' subject ' num2str(s) ' ID:' num2str(this)])
        others=view_vec;
        others(s)=[];
        this_rois=load(['rois/rois_' num2str(this) '.mat']);
        thisTS=this_rois.nodeTS;
        save(['isfc_data/view' num2str(view) '/this_' num2str(s) '.mat'],'thisTS');
        othersTS=0;
        for o=1:length(others)
            load(['rois/rois_' num2str(others(o)) '.mat']);
            othersTS=othersTS+zscore(nodeTS);
        end
        othersTS=zscore(othersTS);
        save(['isfc_data/view' num2str(view) '/others_' num2str(s) '.mat'],'othersTS');
    end
end
end

alldata=[];
for s=1:120
	disp(['loading ' num2str(s)])
	load(['rois/rois_' num2str(s) '.mat']);
	alldata(:,:,s)=nodeTS;
end

mkdir('isc')
for r=1:size(alldata,2)
	slice=squeeze(alldata(:,r,:));
	iscmat=corr(slice);
	save(['isc/iscroi_' num2str(r) '.mat'],'iscmat')
end

function roi_extract_helper(s,subjpath)
    disp(['parfor Rois extraction for subject ' num2str(s)])
    load rois_Power264
    tic
    cfg=[];
    cfg.rois=rois;
    cfg.infile=subjpath{s};
    cfg.write=0;
    cfg.usemean=1;
    % extract rois
    [nodeTS perc]=bramila_roiextract(cfg);
    outfile=['rois/rois_' num2str(s) '.mat'];
    disp(['Storing rois as ' outfile])
    save(outfile,'nodeTS','perc');
    toc
end

