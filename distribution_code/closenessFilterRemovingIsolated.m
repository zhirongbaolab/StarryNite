%discard maxima that are too close to eachother
%use values to discard weaker of them

function [filtpoints,indicies]=closenessFilterRemovingIsolated(points,values,threshold,zscalingfactor,coverage)
scaledpoints=[points(:,1),points(:,2),points(:,3)*zscalingfactor];
valid=zeros(size(values));


sizes=size(points);
%filtpoints=[];
%distances=distance_mem(scaledpoints',scaledpoints'); 
distances=distance(scaledpoints',scaledpoints'); 

for i=1:sizes(1)
   I=find(distances(i,:)<threshold);
   I2=find(I~=i);
   vals2=I(I2);
   %i is largest maxima within distance threshold of itself
   %or is (arbitrarily) the last of a set of equally valued maxima
   if (isempty(I2))
       %its bad if nothing else around
   else
       %its a center if it has good coverage, and is brigther than good
       %coverage around it
       %maxcompetitor=max(values(vals2(coverage(vals2)>13)));
       % if( coverage(i)>13&&(isempty(maxcompetitor)||values(i)>=maxcompetitor))
       if( values(i)>=max(values(vals2)))
    %filtpoints=[filtpoints;points(i,:)];
    valid(i)=1;
       end
   end
end
indicies=logical(valid);
filtpoints=points(indicies,:);


