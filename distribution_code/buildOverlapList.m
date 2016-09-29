function [overlaplist,merges]=buildOverlapList(centers,logodds,range,xycenterslength)

overlaplist=cell(xycenterslength,1);
for i=1:length(centers)
    c=centers{i};
    %l=logodds{i};
    r=range{i};
    for j=1:length(r)
        %add this nucleus to the claimants list for all disks it claims
        overlaplist{c(r(j),4)}=[overlaplist{c(r(j),4)},i];
    end
end
merges=cell(length(centers),1);
for i=1:length(overlaplist)
    ol=overlaplist{i};
    if (length(ol)>1)%more than one nucleus has claim on this disk
        for h=1:length(ol)
            for j=h+1:length(ol)
                one=max(ol(h),ol(j));
                two=min(ol(h),ol(j));
                if(isempty(find(merges{one}==two)))
                    merges{one}=[merges{one},two];
                    %'inside'
                end
            end
        end
    
    end
end