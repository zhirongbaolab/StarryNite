function [sscore,range1,range2]=splitScoreBetweenOverlap(n1,n2,n1c,n2c,centers,ranges,logodds)
%score overlap assuming n1 is the lower numbered planes
c1=centers{n1};
c2=centers{n2};
l1=logodds{n1};
l2=logodds{n2};
rangemask1=zeros(size(l1));
rangemask2=zeros(size(l2));
rangemask1(ranges{n1})=1;
rangemask2(ranges{n2})=1;

%each gets by definition the side closer to them, everything in middle
    %is split in such a way that the half each gets has greatest score
    sscore=sum(find(l1(rangemask1'&c1(:,3)<=n1c(3))));%portion definitely 1 by constraint on topolgoy
    range1=find(rangemask1'&c1(:,3)<=n1c(3));
    sscore=sscore+sum(find(l2(rangemask2'&c2(:,3)>=n2c(3))));%portion definitely 2 by constraint on topolgoy+
    range2=find(rangemask2'&c2(:,3)>=n2c(3));
    
    uncertainrange1=find(rangemask1'&c1(:,3)>n1c(3)&c1(:,3)<n2c(3));
    for i=1:length(uncertainrange1)
        locationin2=find(c1(uncertainrange1(i),4)==c2(:,4).*rangemask2');%find if this disk # appears in others range
        if(~isempty(locationin2))
            if(l2(locationin2)<=l1(uncertainrange1(i)))
                sscore=sscore+l1(uncertainrange1(i));
                range1=[range1;uncertainrange1(i)];
            end
        else
            sscore=sscore+l1(uncertainrange1(i)); %since its not in other range I keep it
            range1=[range1;uncertainrange1(i)];
        end
    end
    uncertainrange2=find(rangemask2'&c2(:,3)<n2c(3)&c2(:,3)>n1c(3));
    for i=1:length(uncertainrange2)
        locationin1=find(c2(uncertainrange2(i),4)==c1(:,4).*rangemask1');%find if this disk # appears in others range
        if(~isempty(locationin1))
            if(l1(locationin1)<l2(uncertainrange2(i)))
                sscore=sscore+l2(uncertainrange2(i));
                 range2=[range2;uncertainrange2(i)];
            end
        else
            sscore=sscore+l2(uncertainrange2(i)); %since its not in other range I keep it
            range2=[range2;uncertainrange2(i)];
        end
    end
    