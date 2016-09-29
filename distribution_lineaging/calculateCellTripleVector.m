function [ datavector ] =calculateCellTripleVector(nucs,i,nucs2,j,nucs3,k,anisotropy,interval)
%{
%score measurements related to agreement between daughters

midpoint=(nucs3.finalpoints(k,:)+nucs2.finalpoints(j,:))./2;
drift=distance_anisotropic(nucs.finalpoints(i,:)',midpoint',anisotropy);
daughterdist=distance_anisotropic(nucs3.finalpoints(k,:)',nucs2.finalpoints(j,:)',anisotropy);
if(daughterdist==0)
'weird';
end
%ensure that j is smaller of distances
%distj=distance_anisotropic(nucs.finalpoints(i,:)',nucs2.finalpoints(j,:)',anisotropy);
%distk=distance_anisotropic(nucs.finalpoints(i,:)',nucs3.finalpoints(k,:)',anisotropy);

%if(distj>distk)
%temp=k;
%tempnucs=nucs3;
%k=j;
%nucs3=nucs2;
%j=temp;
%nucs2=tempnucs;
%end
 % min(distj/distk,distk/distj),...
    
  %     distj/distk,...
 datavector=[log(drift./daughterdist+1),...
    nucs2.totalGFP(j)/nucs3.totalGFP(k),...
      nucs2.avgGFP(j)/nucs3.avgGFP(k),...
    nucs2.finaldiams(j)/nucs3.finaldiams(k)];

end
%}
%calculate new motion cues
%


midpoint=(nucs3.finalpoints(k,:)+nucs2.finalpoints(j,:))./2;
midpoint(3)=midpoint(3).*anisotropy(3);
divline=nucs2.finalpoints(j,:)-nucs3.finalpoints(k,:);
divline(3)=divline(3).*anisotropy(3);
otherline=nucs.finalpoints(i,:)-nucs3.finalpoints(k,:);
otherline(3)=otherline(3).*anisotropy(3);
%proj of parent onto division axis
p1=nucs3.finalpoints(k,:);
p1(3)=p1(3)*anisotropy(3);
parentproj=p1+(dot(otherline,divline)/dot(divline,divline))*divline;

perpdrift=distance(nucs.finalpoints(i,:)',parentproj');
paralleldrift=distance(midpoint',parentproj');
xycomponent=(divline(1)^2+divline(2)^2)^.5;
zcomponent=abs(divline(3));

datavector=[perpdrift,paralleldrift,xycomponent,zcomponent];
normv=mean(nucs.selfdistance);
if(isinf(normv))
    normv=mean(nucs.finaldiams);
end

datavector=datavector./normv;
datavector=datavector./interval;
datavector=[datavector, nucs2.totalGFP(j)/nucs3.totalGFP(k),...
      nucs2.avgGFP(j)/nucs3.avgGFP(k)];

%compute aspect ratio from segmentation

global parameters
slices=nucs.merged_sliceindicies{i};
positions=[];
filterslice=true;
if(filterslice)
    maxvalue=max(nucs.diskintensity(slices));
    slices=slices(nucs.diskintensity(slices)>=maxvalue*parameters.boundary_percent);
end

for s=1:length(slices)
    %nned to add in z still
    % positions=[positions;result.xpositions(:,slices(s)),result.ypositions(:,slices(s)),repmat(result.sliceCenters(slices(s),3),16,1)];
    nonzero=find(nucs.xpositions(:,slices(s))~=0|nucs.ypositions(:,slices(s))~=0);
    positions=[positions;nucs.xpositions(nonzero,slices(s)),nucs.ypositions(nonzero,slices(s)),repmat(nucs.sliceCenters(slices(s),3)*anisotropy(3),length(find(nonzero)),1)];
end
if(length(slices)==1)
    positions=[positions;positions+repmat([0,0,anisotropy(3)],size(positions,1),1)];
end
[coeff,score,roots] = princomp(positions);
[coeff2,score2,roots2] = princomp(positions(:,1:2));

%add aspect ratio of parent to vector (3d, 2d) 
%datavector=[datavector,roots(2)/roots(3),roots2(1)/roots2(2)];
datavector=[datavector,roots2(1)/roots2(2)];

slices=nucs2.merged_sliceindicies{j};
positions=[];
filterslice=true;
if(filterslice)
    maxvalue=max(nucs2.diskintensity(slices));
    slices=slices(nucs2.diskintensity(slices)>=maxvalue*parameters.boundary_percent);
end

for s=1:length(slices)
    %nned to add in z still
    % positions=[positions;result.xpositions(:,slices(s)),result.ypositions(:,slices(s)),repmat(result.sliceCenters(slices(s),3),16,1)];
    nonzero=find(nucs2.xpositions(:,slices(s))~=0|nucs2.ypositions(:,slices(s))~=0);
    positions=[positions;nucs2.xpositions(nonzero,slices(s)),nucs2.ypositions(nonzero,slices(s)),repmat(nucs2.sliceCenters(slices(s),3)*anisotropy(3),length(find(nonzero)),1)];
end
if(length(slices)==1)
    positions=[positions;positions+repmat([0,0,anisotropy(3)],size(positions,1),1)];
end

[coeffd1,scored1,rootsd1] = princomp(positions);
[coeff2d1,score2d1,roots2d1] = princomp(positions(:,1:2));

slices=nucs3.merged_sliceindicies{k};
positions=[];
filterslice=true;
if(filterslice)
    maxvalue=max(nucs3.diskintensity(slices));
    slices=slices(nucs3.diskintensity(slices)>=maxvalue*parameters.boundary_percent);
end

for s=1:length(slices)
    %nned to add in z still
    % positions=[positions;result.xpositions(:,slices(s)),result.ypositions(:,slices(s)),repmat(result.sliceCenters(slices(s),3),16,1)];
    nonzero=find(nucs3.xpositions(:,slices(s))~=0|nucs3.ypositions(:,slices(s))~=0);
    positions=[positions;nucs3.xpositions(nonzero,slices(s)),nucs3.ypositions(nonzero,slices(s)),repmat(nucs3.sliceCenters(slices(s),3)*anisotropy(3),length(find(nonzero)),1)];
end
if(length(slices)==1)
    positions=[positions;positions+repmat([0,0,anisotropy(3)],size(positions,1),1)];
end

[coeffd2,scored2,rootsd2] = princomp(positions);
[coeff2d2,score2d2,roots2d2] = princomp(positions(:,1:2));



% d1, d2 change in 3d d1,d2 change in 2d
%datavector=[datavector,rootsd1(2)/roots(2),rootsd2(2)/roots(2),roots2d1(1)/roots2(1),roots2d2(1)/roots2(1)];

%3d, 2d daughter size agreement
%datavector=[datavector,rootsd1(2)/rootsd2(2),roots2d1(1)/roots2d2(1)];

datavector=[datavector,roots2d1(1)/roots2(1),roots2d2(1)/roots2(1)];

%3d, 2d daughter size agreement
datavector=[datavector,roots2d1(1)/roots2d2(1)];

