%discard maxima that are too close to eachother
%use values to discard weaker of them

function [filtpoints,indicies]=closenessFilter(points,values,threshold,zscalingfactor)
scaledpoints=[points(:,1),points(:,2),points(:,3)*zscalingfactor];
valid=zeros(size(values));


sizes=size(points);
filtpoints=[];
distances=distance(scaledpoints',scaledpoints'); 

for i=1:sizes(1)
   I=find(distances(i,:)<threshold);
   I2=find(I~=i);
   vals2=I(I2);
   %i is largest maxima within distance threshold of itself
   %or is (arbitrarily) the last of a set of equally valued maxima
   %if (length(I2)==0 || values(i)>max(values(vals2)) || (values(i)==max(values(vals2))&i>=max(I)))
   if (length(I2)==0 || values(i)>+max(values(vals2)) )

    %filtpoints=[filtpoints;points(i,:)];
    valid(i)=1;
   end
end
indicies=logical(valid);
filtpoints=points(indicies,:);


