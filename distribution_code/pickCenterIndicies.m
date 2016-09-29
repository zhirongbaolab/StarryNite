function centerindicies=pickCenterIndicies(X,twodmax,celldiameter,anisotropy,numcells)
%3d max

sizes=size(X);

%3d max implemented manually
%offsets into 18 3D neighbors not examined in computing 2d max
offsets=[-1,-1;-1,0;-1,1;0,-1;0,0;0,1;1,-1;1,0;1,1];
offsets=[offsets,ones(length(offsets),1);offsets,-1*ones(length(offsets),1)];
[x,y,z]=ind2sub(sizes,twodmax);
twodxy=[x,y,z];

%expanded2dxyz=repmat(2dxyz,
centerindicies=[];
for i=1:length(twodmax)
    points=[twodxy(i,1)+offsets(:,1),twodxy(i,2)+offsets(:,2),twodxy(i,3)+offsets(:,3)];
    points(:,:)=max(points(:,:),1);
    points(:,1)=min(points(:,1),sizes(1));
    points(:,2)=min(points(:,2),sizes(2));
     points(:,3)=min(points(:,3),sizes(3));
     %now that have subscritp indicies of neighbors turn to linear
     neighbors=sub2ind(sizes,points(:,1),points(:,2),points(:,3));
     if (max(X(neighbors))<=X(twodmax(i)))
         centerindicies=[centerindicies;i];
     end
end


%{
tic
%3d max version
%apply zlevel bottom threshold
X(:,:,zlevel)=X(:,:,zlevel-1.)/2;

X(:,:,zlevel+1:sizes(3))=0;

Xmax=imregionalmax(X,26);
% locations of 3d maxima indicies into image
maxima=find(Xmax);
maxima=maxima(X(maxima)>maximathreshold);

%translate into indicies into xy maximas
centerindicies=[];
for maxi=1:length(maxima)
       centerindicies=[centerindicies;find(twodmax==maxima(maxi))];
end
toc
%}

%{
%filtered xy max version
[altxymax,indicies]=closenessFilter(twodxy,X(twodmax),getParameter('selection_dist',numcells)*celldiameter,anisotropy);
centerindicies=find(indicies);
%}