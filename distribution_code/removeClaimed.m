function [unclaimed,claimed]=removeClaimed(centers,ranges,xymax)
%return unclaimed xy points 
s=size(xymax);
claimed=zeros(1,s(1));

 for i=1:length(centers)


  points=centers{i};
  claimed(points(ranges{i},4))=1;
    
 end
 unclaimed=xymax((~claimed),:);
 claimed=~claimed;
 
 