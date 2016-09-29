function bestrange=vcalculateMaximalRange(planes,centerpoint,logodds,xycoverage, numcells)
%first negative disk cut with fudge factor

badthreshold=13; %coverage threshold for considering a disk unreliable

%start with all center points valid
%their logodds scores are meaningless anyway
range=find(planes(:,3)==centerpoint(3))';
bestrange=range;
%score=0;
%bestscore=0;

pivot=centerpoint(3);
topplane=min(min(planes(:,3)),pivot);
bottomplane=max(max(planes(:,3)),pivot);

%ammount by which a range can be worse than the best range and still be
%used because it is bigger
%threshold=1;
threshold=getParameter('rangethreshold',numcells,centerpoint);

pos=true;
%loop over possible top side sets
for j=pivot-1:-1:topplane
    currentdisks=find(planes(:,3)==j);
    % 'good disks' are those above attachment threshold or with coverage so
    % low we assume them to be part of the range
    positive=find(logodds(currentdisks)>-threshold|xycoverage(planes(currentdisks,4))<=badthreshold);
    % positive=find(logodds(currentdisks)>-threshold);
    if(~isempty(positive)&pos==true) %positive disks here so include and continue
        range=[range,currentdisks(positive)'];
      %  score=score+sum(logodds(currentdisks(positive)));

       % bestscore=score;
        bestrange=range;

    else %done

        pos=false;
    end

end

pos=true;
range=bestrange;%roll back worse possibilities considered after best
%score=bestscore;
%loop over possible bottom side sets
for j=pivot+1:1:bottomplane
    currentdisks=find(planes(:,3)==j);
    %  positive=find(logodds(currentdisks)>-threshold);

    positive=find(logodds(currentdisks)>-threshold|xycoverage(planes(currentdisks,4))<=badthreshold);
    if(~isempty(positive)&pos==true) %positive disks here so include and continue
        range=[range,currentdisks(positive)'];
        %score=score+sum(logodds(currentdisks(positive)));

        %bestscore=score;
        bestrange=range;

    else %done
        pos=false;
    end

end
