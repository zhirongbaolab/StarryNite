%discover the set of x,y maximas that belong (are likely caused by based on our model)
%to the center in center
function centers=assign_planes(center,diam,volsize,xymaximas,anisotropy)
distances=distanceToPoint(center(1:2)',xymaximas(:,1:2)');
%distances=[];
centers_up=assignPlanesDirectional(distances,center,diam,volsize,xymaximas,-1,anisotropy);
centers_down=assignPlanesDirectional(distances,center,diam,volsize,xymaximas,+1,anisotropy);
centers=assignPlanesDirectional(distances,center,diam,volsize,xymaximas,0,anisotropy);

centers=[centers_up;centers;centers_down];
