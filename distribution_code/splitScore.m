function [sscore,range1,range2]=splitScore(n1,n2,n1c,n2c,centers,ranges,logodds)
c1=centers{n1};
c2=centers{n2};
l1=logodds{n1};
l2=logodds{n2};
rangemask1=zeros(size(l1));
rangemask2=zeros(size(l2));
rangemask1(ranges{n1})=1;
rangemask2(ranges{n2})=1;

%special case 2 share center then each gets whichever half of rante
%maximizes total score
if(n1c(3)==n2c(3))
    sscore=sum(l1((rangemask1'&c1(:,3)>n2c(3))))+sum(l2((rangemask2'&c2(:,3)<n2c(3))));
    range1=find(rangemask1'&c1(:,3)>=n2c(3));
    range2=find(rangemask2'&c2(:,3)<=n2c(3));
    
    sscore2=sum(l2((rangemask2'&c2(:,3)>n2c(3))))+sum(l1((rangemask1'&c1(:,3)<n2c(3))));
    if(sscore2>sscore)
            range1=find(rangemask1'&c1(:,3)<=n2c(3));
            range2=find(rangemask2'&c2(:,3)>=n2c(3));     
    end
    sscore=max(sscore,sscore2);
end
if(n1c(3)<n2c(3))

    [sscore,range1,range2]=splitScoreBetweenOverlap(n1,n2,n1c,n2c,centers,ranges,logodds);
  
end

if(n1c(3)>n2c(3))
        %note that because input order is changed, output order is changed
        [sscore,range2,range1]=splitScoreBetweenOverlap(n2,n1,n2c,n1c,centers,ranges,logodds);
end

% average of individual scores rather than split
%sscore=sum([l1(ranges{n1}),l2(ranges{n2})])/2;