function [averagedpoints,points,aspectratios,mdiameters,mmaximas, merged_sliceindicies,mergedlogoddssum]=mergeAndJudgeNuclei(diskSet,nucleiSet,mergesets,anisotropy)
%(centers,ranges,centerpoints,diameters,mergesets,anisotropy,validnuc)
%go through joint 1st and 2nd round list, if merged replace conn of points
%with a single point at average
%if not merged keep point as long as it is not unexpectedly small <2/3 size
%in z given xy diam and xy diam >=2/3 global expected xy diam
points=[];
averagedpoints=[];
aspectratios=[];
mdiameters=[];
mmaximas=[];
merging=[];
mergedlogoddssum=[];
merged_sliceindicies={}; % centers{i}(range{i}) indicies for final nuc components

logodds=nucleiSet.logodds;
centers=nucleiSet.centers;
ranges=nucleiSet.range;
diameters=diskSet.xydetdiameters(nucleiSet.centerindicies);
maxvals=diskSet.xymaximavals(nucleiSet.centerindicies);
centerpoints=diskSet.centeredxymax(nucleiSet.centerindicies,:);



%add merged points to set and keep track of all merged points
for i=1:length(mergesets)
    logoddsset=[];
    curr=mergesets{i};
    merging=[merging;curr];%keep track of everything that has been merged
    pointset=[];
    for j=1:length(curr)
        logoddscurr=logodds{curr(j)};
        cc=centers{curr(j)};
        cc=cc(ranges{curr(j)},:);
        logoddsset=[logoddsset,logoddscurr(ranges{curr(j)})];
        pointset=[pointset;cc];
    end
    points=[points;mean(pointset(:,1:3))];
     averagedpoints=[averagedpoints;mean(pointset(:,1:3))];
    height=max(pointset(:,3))-min(pointset(:,3))+1;
    merged_sliceindicies=[merged_sliceindicies,pointset(:,4)];
    %merged_sliceindicies={merged_sliceindicies{:},pointset(:,4)};
    aspectratios=[aspectratios;height*anisotropy/max(diameters(curr))];
    mdiameters=[mdiameters;max(diameters(curr))];
    mmaximas=[mmaximas;max(maxvals(curr))];
    mergedlogoddssum=[mergedlogoddssum;sum(logoddsset)];
end

%currently not discarding anything
for i=1:length(centers)
    if(isempty(find(i==merging, 1)))%if not merged add to point set
        cc=centers{i};
        cc=cc(ranges{i},:);
        height=max(cc(:,3))-min(cc(:,3))+1;
        logoddscurr=logodds{i};
        logoddsset=logoddscurr(ranges{i});
        points=[points;centerpoints(i,:)];
        if(size(cc,1)==1)
            averagedpoints=[averagedpoints;cc(:,1:3)];
        else
            averagedpoints=[averagedpoints;mean(cc(:,1:3))];
        end
        aspectratios=[aspectratios;height*anisotropy/diameters(i)];
        mdiameters=[mdiameters;diameters(i)];
        mmaximas=[mmaximas;maxvals(i)];
        merged_sliceindicies=[merged_sliceindicies ,cc(:,4)];
       %merged_sliceindicies={merged_sliceindicies{:},cc(:,4)};
        mergedlogoddssum=[mergedlogoddssum;sum(logoddsset)];
    end
end