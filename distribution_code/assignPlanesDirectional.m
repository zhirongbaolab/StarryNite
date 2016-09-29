
%discover the set of x,y maximas in cylinder around center
function centers=assignPlanesDirectional(distances,center,diameter,Xfilt,xymaximas,zoffset,xydetdiameters,anisotropy)
volsize=size(Xfilt);
centers=[];
  querydiameter=diameter;

potential=find(xymaximas(:,3)==center(3)+zoffset);
%potential=xymaximas(:,3)==(center(3)+zoffset);
%potential=potential&(xymaximas(:,1)<(center(1)+diameter/2));
%potential=potential&(xymaximas(:,1)>(center(1)-diameter/2));
%potential=potential&(xymaximas(:,2)<(center(2)+diameter/2));
%potential=potential&(xymaximas(:,2)>(center(2)-diameter/2));
%potential=find(potential);

%%dist=distance(xymaximas(potential,:)',[center(1),center(2),center(3)+zoffset]');
%dist=distanceToPoint([center(1),center(2),center(3)+zoffset]',xymaximas(potential,:)');
dist=distances(potential);
smalldist=(dist<querydiameter/2);
possiblemaximas=xymaximas(potential(smalldist),:);
%possiblemaximas=xymaximas(potential,:);
%possiblemaximas=possiblemaximas(smalldist,:);

%indicies=potential(smalldist);
indicies=potential(smalldist);

sizes=size(possiblemaximas);
for i=1:sizes(1)
    querycenter=possiblemaximas(i,:);

   if (abs(zoffset*anisotropy/(diameter/2))<3) %short circut to return whole collumn
       
               centers=[centers;[querycenter,indicies(i)]];
    end
end
 nextplane=center(3)+zoffset+sign(zoffset);
 
if(zoffset~=0&&~isempty(centers)&&nextplane>=1&&nextplane<=volsize(3)) %if found something here and next plane is in volume
       centers=[centers;assignPlanesDirectional(distances,center,diameter,Xfilt,xymaximas,zoffset+sign(zoffset),xydetdiameters,anisotropy)];
end
    