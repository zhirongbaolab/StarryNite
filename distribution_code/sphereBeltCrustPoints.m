%indices of all sampled points insides anisotropic in z sphere  
function indicies=sphereBeltCrustPoints(cellsize,centroid,anisotropy,Xorig)

        R = cellsize/2.0;
        h = ones(round(cellsize),round(cellsize),round(cellsize/anisotropy));
        siz=[cellsize+2,cellsize+2,round(cellsize/anisotropy)+2];
        siz = round((siz-1)/2);
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        I = (x.*x/R^2+y.*y/R^2+z.*z/(R/anisotropy)^2);
        
        [f,v]=isosurface(I,1);
        v=round(v);
        
        x=v(:,2);
        y=v(:,1);
        z=v(:,3);
        x=x-siz(2)-1+centroid(1);
        y=y-siz(1)-1+centroid(2);
        z=z-siz(3)-1+centroid(3);
        
        %short circut to return only center plane
        centerI=find(z==centroid(3));
        z=z(centerI);
        x=x(centerI);
        y=y(centerI);
        
        
        sizes=size(Xorig);
           z=min(z,sizes(3)); %clamping hack  if spheres wind up outside volume doesnt really matter what do since these are never cells anyway
          z=max(z,1);
        y=max(y,1);
        y=min(y,sizes(1));
        x=max(x,1);
        x=min(x,sizes(2));
        
        indicies=sub2ind(size(Xorig),y,x,z);