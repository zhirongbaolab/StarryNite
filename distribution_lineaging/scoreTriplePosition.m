function [ minscore,configuration] = scoreTriplePosition( points1,points2,anisotropy )
% test each possible pairing of the 3 points in points1 and points 2 via
% summed distance,return distance and best configuration
points1(:,3)=points1(:,3)*anisotropy;
points2(:,3)=points2(:,3)*anisotropy;

scores=[];
configurations=[];
for i=1:3
    for j=1:3
        for k=1:3
            if(i~=j&i~=k&j~=k)
                configurations=[configurations;i,j,k];
                %sum of distances for all 3 pairings
                score=distance(points1(1,:)',points2(i,:)')+...
                    distance(points1(2,:)',points2(j,:)')+...
                    distance(points1(3,:)',points2(k,:)');
                scores=[scores;score];
            end
        end
    end
end
[minscore,i]=min(scores);
configuration=configurations(i,:);

end

