function mscore=mergeScore2(n1,n2,n1c,n2c,centers,ranges,xymax,xymaximavals,xydetdiameters,anisotropy,xycoverage)
%pull out unique disks and score their union for logodds
c1=centers{n1};
c2=centers{n2};
mergedset=[c1(ranges{n1},:);c2(ranges{n2},:)];
[u,i,j]=unique(mergedset(:,4));
mergedset=mergedset(i,:);
%set center to geometric center of union range
centerplane=round(mean([n1c(3),n2c(3)]));
%take xy maxima there as new center if more than one take first
centerindex=mergedset(find(mergedset(:,3)==centerplane,1),4);
%arbitrarily take the first if there happens to be two on that plane
%if(length(centerindex>1))
%    centerindex=centerindex(1);
%end
%'remember not using prediction in logodds calculation for merging part'
mlogodds=calculateLogodds(mergedset,xymax(centerindex,:),xymaximavals(centerindex),xydetdiameters(centerindex),xymaximavals,xydetdiameters,anisotropy,xycoverage);
mscore=sum(mlogodds);
