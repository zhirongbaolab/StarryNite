function diskSet=calculateSliceGFP(Xorig,diskSet)
global parameters;
diskSet.GFPsums=zeros(1,length(diskSet.xydetdiameters));
diskSet.diskArea=zeros(1,length(diskSet.xydetdiameters));
diskSet.diskMax=zeros(1,length(diskSet.xydetdiameters));
raysort2s=[1,13,5,9,3,12,7,15,2,14,6,10,4,11,8,16];
for i=1:length(diskSet.xydetdiameters)
    
    
    
    good=logical(diskSet.xpositions(:,i)~=0|diskSet.ypositions(:,i)~=0);
    
    sortedgood=good(raysort2s);
    sortedx=diskSet.xpositions(raysort2s,i);
    sortedy=diskSet.ypositions(raysort2s,i);
     sx=size(Xorig);
 
    
    
    %need to order and filter points
    %
    %old
  %  pointsmask=poly2mask(sortedx,sortedy,size(Xorig,1),size(Xorig,2));
  %  points=find(pointsmask);
     % [x1,x2]=ind2sub(sx(1:2),points);
     
     %   linind=sub2ind(sx,x1,x2,ones(size(x1)).*diskSet.xymax(i,3));
    %new
    minx=floor(min(sortedx));
    maxx=ceil(max(sortedx));
    miny=floor(min(sortedy));
    maxy=ceil(max(sortedy));
    testx=repmat(linspace(minx,maxx,maxx-minx+1),1,maxy-miny+1); 
    testy=repmat(linspace(miny,maxy,maxy-miny+1),maxx-minx+1,1);
    testy=reshape(testy,1,numel(testy));
    in=inpolygon(testx,testy,sortedx,sortedy);
    inpsmall=find(in);
    [x1s,x2s]=ind2sub([maxx-minx+1,maxy-miny+1],inpsmall);
    x1s=x1s+min(testx)-1; %not entirely sure why there is a one off error here but quite clearly is
    x2s=x2s+min(testy)-1;
    %x1s should match x1 if this is right
    %passed in in opposite order to old because computed in geometric not
    %raster space
    linind=sub2ind(sx,x2s,x1s,ones(size(x2s)).*diskSet.xymax(i,3));
   
        %{
    diskSet.diskArea(i)=length(points);
    slice=Xorig(:,:,diskSet.xymax(i,3));
    diskSet.GFPsums(i)=sum(slice(points));
    %}
    %attempt a more efficient version of above which is hotspot in code
    %diskSet.GFPsums(i)=sum(sum(pointsmask.*Xorig(:,:,diskSet.xymax(i,3))));
     diskSet.GFPsums(i)=sum(Xorig(linind));
    
    
    %{
    %compute median value of non background
    %I cant find where this is used so I'm trying to rem it out
    if(~isfield(parameters,'polarBackgroundThreshold'))
       % points=find(pointsmask&slice>2200);
       points=find(pointsmask);
    else
        points=find(pointsmask&Xorig(:,:,diskSet.xymax(i,3))>parameters.polarBackgroundThreshold);
    end
    
    if(isempty(points))
        diskSet.diskMax(i)=nan;
    else
    diskSet.diskMax(i)=prctile(slice(points),95);
    end
    %}
    
  %compute median value of non background

    if(~isfield(parameters,'polarBackgroundThreshold'))
        %  points=find(pointsmask);
        %dont need to do anything points already good
    else
       % points=find(pointsmask&Xorig(:,:,diskSet.xymax(i,3))>parameters.polarBackgroundThreshold); 
       
       %[x1,x2]=ind2sub(sx(1:2),points);
       % linind=sub2ind(sx,x1,x2,ones(size(x1)).*diskSet.xymax(i,3));
       %filter without having explicit mask
       linind=linind(Xorig(x1,x2,diskSet.xymax(i,3))>parameters.polarBackgroundThreshold);
    end
    
    if(isempty(linind))
        diskSet.diskMax(i)=nan;
    else
        
        %diskSet.diskMax(i)=prctile(slice(points),95);
    diskSet.diskMax(i)=prctile(Xorig(linind),95);
    
    end
    
    
end
diskSet;