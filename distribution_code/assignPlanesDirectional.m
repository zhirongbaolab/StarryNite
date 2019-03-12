
%discover the set of x,y maximas in cylinder around center
function centers=assignPlanesDirectional(distances,center,diameter,volsize,xymaximas,zoffset,anisotropy)
%volsize=size(Xfilt);
centers=[];
  querydiameter=diameter;

%this block is takinig 1/10th of execution time

potential=find(xymaximas(:,3)==center(3)+zoffset);
dist=distances(potential);
smalldist=(dist<querydiameter/2);
possiblemaximas=xymaximas(potential(smalldist),:);
indicies=potential(smalldist);

%{
%is this faster?
%no%
possiblieindicies=(xymaximas(:,3)==center(3)+zoffset&distances'<querydiameter/2);
possiblemaximas=xymaximas(possiblieindicies,:);
indicies=find(possiblieindicies);
%}
sizes=size(possiblemaximas);
for i=1:sizes(1)
    querycenter=possiblemaximas(i,:);

   if (abs(zoffset*anisotropy/(diameter/2))<3) %short circut to return whole collumn
       
               centers=[centers;[querycenter,indicies(i)]];
    end
end
 nextplane=center(3)+zoffset+sign(zoffset);
 
if(zoffset~=0&&~isempty(centers)&&nextplane>=1&&nextplane<=volsize(3)) %if found something here and next plane is in volume
       centers=[centers;assignPlanesDirectional(distances,center,diameter,volsize,xymaximas,zoffset+sign(zoffset),anisotropy)];
end
    