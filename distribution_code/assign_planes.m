%discover the set of x,y maximas that belong (are likely caused by based on our model)
%to the center in center
function centers=assign_planes(center,diam,Xfilt,xymaximas,xydetdiameters,anisotropy)
distances=distanceToPoint(center(1:2)',xymaximas(:,1:2)');
%distances=[];
centers_up=assignPlanesDirectional(distances,center,diam,Xfilt,xymaximas,-1,xydetdiameters,anisotropy);
centers_down=assignPlanesDirectional(distances,center,diam,Xfilt,xymaximas,+1,xydetdiameters,anisotropy);
centers=assignPlanesDirectional(distances,center,diam,Xfilt,xymaximas,0,xydetdiameters,anisotropy);

centers=[centers_up;centers;centers_down];
