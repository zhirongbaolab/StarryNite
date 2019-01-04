function [trackingparameters,esequence]=initializeTrackingStructures(esequence,trackingparameters)
%initialize data structures and perform precomputations for tracking 

%bring global detection parameters in to allow gfp computations
global parameters;

%initialize structures for sucessors and deletion
for t=trackingparameters.starttime:trackingparameters.endtime
    if(isempty(esequence{t}))
        esequence{t}.finalpoints=[];
    end
    esequence{t}.suc=-1*ones(size(esequence{t}.finalpoints,1),2);
    esequence{t}.suc_time=-1*ones(size(esequence{t}.finalpoints,1),2);
    esequence{t}.pred=-1*ones(size(esequence{t}.finalpoints,1),1);
    esequence{t}.pred_time=-1*ones(size(esequence{t}.finalpoints,1),1);
    
    esequence{t}.delete=zeros(size(esequence{t}.finalpoints,1),1);
end


for i=trackingparameters.starttime:trackingparameters.endtime
    if(~isempty(esequence{i}.finalpoints))
    %initialize avg nn distance
    distances=distance_anisotropic(esequence{i}.finalpoints',esequence{i}.finalpoints',trackingparameters.anisotropyvector);
    for j=1:max(size(distances))
        distances(j,j)=Inf;
    end
    mindistances=min(distances);
    diam=mean(mindistances);
    
     trackingparameters.forwardcutoff(i)=diam*trackingparameters.candidateCutoff;%1.5;%candidate cutoff (distance)
     esequence{i}.selfdistance=mindistances;
    
    
    %integrate GFP
    [integratedGFP,area,maxslice,diskMax]=integrateGFP(esequence{i},parameters);
    esequence{i}.totalGFP=integratedGFP;
    esequence{i}.avgGFP=integratedGFP./area;
    esequence{i}.maxslice=maxslice;
    esequence{i}.mdiskMax=diskMax;
    else
        trackingparameters.forwardcutoff(i)=-1;%1.5;%candidate cutoff (distance)
        esequence{i}.selfdistance=[];
        esequence{i}.totalGFP=[];
        esequence{i}.avgGFP=[];
    end
end


if (isfield(trackingparameters,'abscutoff')&&trackingparameters.abscutoff)
    'absolute tracking working working'
    trackingparameters.forwardcutoff=trackingparameters.candidateCutoff*ones(size(trackingparameters.forwardcutoff));
end


%calculate confidences
for time=trackingparameters.starttime:trackingparameters.endtime
       if(~isempty(esequence{time}.finalpoints))

    confidencedata=calculateConfidenceVector(esequence{time},parameters);
    esequence{time}.confidencevector=confidencedata;
    tempconfidences=mvnpdf(confidencedata, trackingparameters.model.con_goodmean,trackingparameters.model.con_goodstd)...
        ./mvnpdf(confidencedata, trackingparameters.model.con_badmean,trackingparameters.model.con_badstd);
    tempconfidences=1./tempconfidences; %pdf ratio of bad
    %use raw confidences
    tempconfidences(isnan(tempconfidences))=Inf;%...
    esequence{time}.confidences=tempconfidences;
       else
            esequence{time}.confidencevector=[];
                esequence{time}.confidences=[];
       end
       
end

