function mergegroups=conFromPairMerge(mergeinfo)
%given list of pairwise merges in rows unravel into cell array of connected
%components of all n things merged
mergegroups={};
s=size(mergeinfo);
counter=0;
for i=1:s(1)

    m1=mergeinfo(i,1);
    m2=mergeinfo(i,2);
    if(m1~=-1)
        merges=[];
        %indicies of everyplace that m1 or m2 occur
        rows=find(mergeinfo(:,1)==m1|mergeinfo(:,2)==m1|mergeinfo(:,2)==m2|mergeinfo(:,1)==m2);
        while(~isempty(rows))
            merges=[merges;unique([mergeinfo(rows,1);mergeinfo(rows,2)])];
            mergeinfo(rows,:)=-1; %wipe out merges already procesed
            rows=zeros(s(1),1);
            for j=1:length(merges)
                rows=rows|mergeinfo(:,1)==merges(j)|mergeinfo(:,2)==merges(j);
            end
            rows=find(rows);
        end
            
    mergegroups{counter+1}=merges;
    counter=counter+1;
    end

end
