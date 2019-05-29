function [diskSet,centerindicies]=createDiskSet(X,maximathreshold,zlevel,celldiameter,anisotropy,numcells,ROIpoints,ROI,ROIxmin,ROIymin)

s=size(X);
diskSet=struct;

%2d max
Xmax2d=zeros(size(X),'int8');
parfor i=1:s(3)
    Xmax2d(:,:,i)=imregionalmax(X(:,:,i));
end

%wipe out bottom XYmaxes before they are even considered anywhere
%Xmax2d(:,:,zlevel:s(3))=0;
twodmax=find(Xmax2d);


twodmax=twodmax(X(twodmax)>maximathreshold);%apply threshold to maximas
%{
%'warning fish hack active'
Xr=imresize(Xr,.5);
for i=1:size(Xr,3)
Xr(:,:,i)=medfilt2(Xr(:,:,i),[13,13],'symmetric');
end
Xr=imresize(Xr,2);

twodmax=twodmax(Xr(twodmax)>1750);%was 1740 more aggressive filtering should allow to rais it
%}
[x,y,z]=ind2sub(s,twodmax);

%filter to ROI region if it exists
if(~isempty(ROIpoints))
        if(~isempty(twodmax))
            pos=[y,x,z];
            polygon_points=ROIpoints;
            if(ROI)
                goodpoints=inpolygon(pos(:,1)+ROIxmin,pos(:,2)+ROIymin,polygon_points(:,1),polygon_points(:,2));
            else
                goodpoints=inpolygon(pos(:,1),pos(:,2),polygon_points(:,1),polygon_points(:,2));
           
            end
            x=x(goodpoints);
            y=y(goodpoints);
            z=z(goodpoints);
            twodmax=twodmax(goodpoints);
        end
end


%xymax=[y,x,z];
diskSet.xymax=[y,x,z];
[diskSet.xydetdiameters,diskSet.centeredxymax,diskSet.xycoverage,xpositions,ypositions]=calculateSphereDiameters_geometric(diskSet.xymax,X,celldiameter,true);
diskSet.xymaximavals=X(twodmax);
diskSet.xpositions=xpositions;
diskSet.ypositions=ypositions;
centerindicies=pickCenterIndicies(X,twodmax,celldiameter,anisotropy,numcells);