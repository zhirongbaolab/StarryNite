%not working**


%offline matching to be run before training_code.m if matching has not been run
%online
%assumes just one embryo in dataset which is typically true of recent runs
for time=startime:endtime
    
e=esequence{time};
    %points=mergedNucleiSet.positions;
    points=e.finalpoints;
    
    mergeinfo2a=e.mergeinfo;
pairmerges=e.pairmerges;

    %read labeled data
[celldata,cellnames]=readnuclei([nucleidir,embryonumber_ed,'\t',num2str(time,'%03d'),'-nuclei']);
%read undedited sn data
%[Ucelldata,Ucellnames]=readnuclei([nucleidir,embryonumber,'_unedited\t',num2str(time,'%03d'),'-nuclei']);

    threshold=.75*e.celldiameter; %for identity of detected and labeled nuclei

    celllocations=celldata(:,4:6);
%Ucelllocations=Ucelldata(:,4:6);

%stats from initial 1st round matching matching
celllocations(:,1:2)=celllocations(:,1:2).*downsample;
%Ucelllocations(:,1:2)=Ucelllocations(:,1:2).*downsample;
%remember arguments to sub2ind are y,x,z
diameters=celldata(:,7)*downsample;
filteredpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies(1:e.firstroundlength),:);


[matches,matchessr]=compareDetectionWRadius_3(filteredpointsswitched,celllocations,celldata(:,7)*downsample*.5,1.5,anisotropy);
FP=filteredpointsswitched((matches==-1),:);
FN=celllocations((matchessr==-1),:);
%pmatches=filteredpointsswitched((matches~=-1),:);

totalcorrectnuclei=length(matchessr);
initialFN=length(FN);
initialFP=length(FP);

e.initialFN=FN;
e.initialFP=FP;

%stats from all 1st and 2nd round points matching
%match all centers
[matches2,matchessr2]=compareDetectionWRadius_3(diskSet.centeredxymax(nucleiSet.centerindicies,:),celllocations,celldata(:,7)*downsample*.5,1.5,anisotropy);

if(~isempty(mergeinfo2a))
    if(time~=firsttimestep)
        allvalid=[allvalid,(matches2(mergeinfo2a(:,1))==-1)|(matches2(mergeinfo2a(:,2))==-1)];
    end
end
e.fmatches=matches2;
e.fmatchesr=matchessr2;
e.all2ndroundFN=celllocations((matchessr2==-1),:);

nrl=[];
for i=1:length(nucleiSet.range)
    nrl=[nrl,length(nucleiSet.range{i})];
end
allpoints=diskSet.centeredxymax(nucleiSet.centerindicies,:);
e.zerodiskcorrect=allpoints((nrl<2&matches2==-1),:);
e.zerodiskwrong=allpoints((nrl<2&matches2>0),:);

e.matches2=matches2;
if (~isempty(pairmerges))
    introducedmergeFN=length(find(matches2(pairmerges(:,1))>0&matches2(pairmerges(:,2))>0));
    e.introducedmergeFN=pairmerges((matches2(pairmerges(:,1))>0&matches2(pairmerges(:,2))>0),:);


    correctmerges=length(find(matches2(pairmerges(:,1))==-1|matches2(pairmerges(:,2))==-1));
    e.correctmerges=pairmerges((matches2(pairmerges(:,1))==-1|matches2(pairmerges(:,2))==-1),:);
else
    e.introducedmergeFN=[];
    e.correctmerges=[];
end
%unmergedFP %could explicitly pull this out but it is implicit in remaining
%FP modulo matching error
%'missed FN' %ditto for this
psizes=size(filteredpointsswitched);
secondroundFP=length(find(find(matches2==-1)>psizes(1)));

e.secondroundFP=allpoints((matches2==-1),:);%all nonmatches
e.secondroundFP=e.secondroundFP((find(matches2==-1)>psizes(1)),:);%just the late ones
e.secondroundFPind=find(matches2==-1);%all nonmatches
e.secondroundFPind=e.secondroundFPind((find(matches2==-1)>psizes(1)));%just the late ones


correct2ndroundnuclei=length(find(find(matches2>0)>psizes(1)));
e.correct2ndroundnuclei=allpoints((matches2>0),:);%all
e.correct2ndroundnuclei=e.correct2ndroundnuclei((find(matches2>0)>psizes(1)),:);%all 2nd round
e.correct2ndroundnucleiind=find(matches2>0);%all
e.correct2ndroundnucleiind=e.correct2ndroundnucleiind((find(matches2>0)>psizes(1)));%all 2nd round

%match final merged and judged points
[matches2,matchessr2]=compareDetectionWRadius_3(points,celllocations,celldata(:,7)*downsample*.5,1.5,anisotropy);
points=round(points);
FP=points((matches2==-1),:);
FN=celllocations((matchessr2==-1),:);
%pmatches=points(find(matches2~=-1),:);
s=size(FN);
finalFN=s(1);
s=size(FP);
finalFP=s(1);
e.finalFN=FN;
e.finalFP=FP;
e.nuclei=length(matchessr2);

end
end