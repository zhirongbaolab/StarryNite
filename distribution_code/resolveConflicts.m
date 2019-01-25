function [mergedNucleiSet,e]=resolveConflicts(nucleiSet,diskSet,e,anisotropy,numcells,celldiameter,savedata)

[overlaplist2a,pmerges2a]=buildOverlapList(nucleiSet.centers,nucleiSet.logodds,nucleiSet.range,length(diskSet.xymaximavals));
[mergeinfo2a,newrange2a]=filter_boundary(nucleiSet.centers,diskSet.centeredxymax(nucleiSet.centerindicies,:),nucleiSet.range,nucleiSet.logodds,pmerges2a,diskSet.centeredxymax,diskSet.xymaximavals,diskSet.xydetdiameters,anisotropy,diskSet.xycoverage);

%calculate ar of union
%compare these range sizes to prediction
%if normal merge or if this merge becomes pairmerge

allr=newrange2a;
alld=diskSet.xydetdiameters(nucleiSet.centerindicies) ;

%calculate AR of actual nuclei even if no prediction need to do it
currentar=zeros(length(alld),1,'single');
for i=1:length(alld);
    centerpoints=nucleiSet.centers{i};
    centerpoints=centerpoints(allr{i},3);
    height=max(centerpoints)-min(centerpoints)+1;
    currentar(i)=height*anisotropy/alld(i);
end


s=size(mergeinfo2a);
mergear=zeros(1,s(1),'single');
for i=1:s(1)
    centerpoints1=nucleiSet.centers{mergeinfo2a(i,1)};
    centerpoints2=nucleiSet.centers{mergeinfo2a(i,2)};
    centerpointsm=[centerpoints1(allr{mergeinfo2a(i,1)},3);centerpoints2(allr{mergeinfo2a(i,2)},3)];
    height=max(centerpointsm)-min(centerpointsm)+1;
    mergear(i)=height*anisotropy/max([alld(mergeinfo2a(i,1)),alld(mergeinfo2a(i,2))]);
end

allpoints=diskSet.centeredxymax(nucleiSet.centerindicies,:);

    %190 dist 7
if(~isempty(mergeinfo2a))
    s=size(mergeinfo2a);
    mergelower=getParameter('mergelower',numcells);
    mergesplit=getParameter('mergesplit',numcells);

    splitthreshold=zeros(s(1),1);
    distmergethreshold=zeros(s(1),1);
    arthreshold=zeros(s(1),1);      
    
    %build arrays of threshold parameters per point based on region location
    for possiblemerge=1:s(1)
        splitthreshold(possiblemerge)=getParameter('split',numcells,allpoints(mergeinfo2a(possiblemerge,1),:));
        distmergethreshold(possiblemerge)=getParameter('nndist_merge',numcells,allpoints(mergeinfo2a(possiblemerge,1),:));
        arthreshold(possiblemerge)=getParameter('armerge',numcells,allpoints(mergeinfo2a(possiblemerge,1),:));
    end
    
    goodmerges=(mergeinfo2a(:,3)<splitthreshold&mergeinfo2a(:,4)>mergelower);
    goodmerges=goodmerges|(mergeinfo2a(:,4)>mergeinfo2a(:,3)*mergesplit);%classic merge split rule
   
    spoints=[allpoints(:,1),allpoints(:,2),allpoints(:,3)*anisotropy];
    distancevals=distance_mem(spoints',spoints');
    distancevals=distancevals./celldiameter; %norm by expected diam
    
    onndist=zeros(s(1),1);
    for j=1:s(1)
        onndist(j)=distancevals(mergeinfo2a(j,1),mergeinfo2a(j,2));
    end
    
    goodmerges=goodmerges|(onndist<distmergethreshold);
    goodmerges=goodmerges|(mergear'<arthreshold);
    
    if(savedata)
    e.goodmerges=goodmerges;
    end
    pairmerges=mergeinfo2a((goodmerges),1:4);
else
    pairmerges=[];
end


mergegroups=conFromPairMerge(pairmerges);
%create final point
%[points,mergedaspectratios,mergeddiams]=mergeAndJudgeNuclei(nucleiSet.centers,newrange2a,diskSet.centeredxymax(nucleiSet.centerindicies,:),diskSet.xydetdiameters(nucleiSet.centerindicies),mergegroups,anisotropy,ones(size(newrange2a)));
[averagedpoints,points,mergedaspectratios,mergeddiams,mergemaximas,merged_sliceindicies,mergedlogoddssum]=mergeAndJudgeNuclei(diskSet,nucleiSet,mergegroups,anisotropy);

if (savedata)
  
e.pairmerges=pairmerges;
e.mergear=mergear;
e.allnewrange=newrange2a;
e.mergeinfo=mergeinfo2a;
e.mergegroups=mergegroups;
e.allaspectratio=currentar;
end
e.mergedlogoddssum=mergedlogoddssum;
mergedNucleiSet=struct;
mergedNucleiSet.averagedpositions=averagedpoints;
mergedNucleiSet.positions=points;
mergedNucleiSet.diameters=mergeddiams;
mergedNucleiSet.aspectratios=mergedaspectratios;
mergedNucleiSet.maximas=mergemaximas;
mergedNucleiSet.merged_sliceindicies=merged_sliceindicies;
