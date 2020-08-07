%read previous labeled timestep to simulate using previous step nuclear
%diameter to set filter size

if(~exist('usestaticdiameter'))
    usestaticdiameter=false;
end

if (nodatause||nodata)

    if(isstruct(previous)&~usestaticdiameter)
        if(~isempty(previous.diams)&&length(previous.diams)>10) %update only if found cells last time otherwise dont do anything
            celldiameter=median(previous.diams);%max(mean(previous.diams),celldiameter*.9);%limit to drop of 10% per TP
        end
    else
        celldiameter=firsttimestepdiam*downsample;
    end
else
    nuclei=[nucleibase,'t',num2str(time-1,'%03d'),'-nuclei'];
    [celldata,cellnames]=readnuclei(nuclei);
    celldiameter=mean(celldata(:,7))*downsample; % size
end

if(isstruct(previous))
    numcells=length(previous.finalpoints);
else
    numcells=firsttimestepnumcells;
end

%replace distribution for 10thround
%if(numcells>351)
%    load(distribution_file2);
%end


sigma=celldiameter*getParameter('sigma',numcells);
%hack for drosophila
%maximathreshold=getParameter('intensitythreshold',time);
maximathreshold=getParameter('intensitythreshold',numcells);




Xorig=X;
sizes=size(X);
%{
fish=true;
if(fish)
        thresimage=Xr;
threshhold=1735; 1755; %1780; %1755;%2600;% 1755*1.3;%eh? 1780;%1780;
thresimage(thresimage>threshhold)=threshhold;
%thresimage=thresimage./1780;
thresimage=thresimage-(min(min(min(thresimage))));
thresimage=thresimage./(max(max(max(thresimage))));
%thresimage=log(thresimage);
thresimage=thresimage.^5;%originally 4
X=thresimage.*X;
%figure;imagesc(max(thresimage,[],3));
end
%}

%----------------------------------------------------------------------
%original code: takes less memory but is 50 times slower in a single thread
%CPU. 85% of the execution is spent here because it uses a for loop over
%chunks and 3D convolution (not separation between Gaussian kernels)

%use tiled if saving memory or by default, otherwise use seprable memory
%intensive all at once call

if (exist('conservememory','var')&&conservememory)
    X=tiledogfilter(X,sigma,anisotropy);
else
%Fernando Amat: this code takes more memory (about 2 times) but is much faster
siz2A=[floor(sigma/2)*2+1,floor(sigma/2)*2+1,max(3,floor(sigma/anisotropy/2)*2+1)]-1;
sigA = siz2A/(4*sqrt(2*log(2)));
siz2B=[floor(sigma*1.6/2)*2+1,floor(sigma*1.6/2)*2+1,max(3,floor(sigma*1.6/anisotropy/2)*2+1)]-1;
sigB = siz2B/(4*sqrt(2*log(2)));

X=imgaussianAnisotropy(single(X),sigA,siz2A)-imgaussianAnisotropy(single(X),sigB,siz2B);
%-----------------------------------------------------------------


end
%{
fish=true;
if (fish)
X(Xr<1755)=0;
end
%}
%figure;imagesc(max(X,[],3));

%if(isnan(max(max(max(X)))))
%    X=zeros(size(X));
%end

if(nodata||nodatause&singlevolume)
    s=size(X);
    for p=1:s(3)
        if(max(max(X(:,:,p)))>maximathreshold*1.5)
            zlevel=p;
        end
    end
else
    %zlevel already defined if data or full sequence(to full volume)
end

if(~exist('ROIpoints'))
    ROIpoints=[];
end

    %disk representation
nucleiSet=struct;
[diskSet,nucleiSet.centerindicies]=...
    createDiskSet(X,maximathreshold,zlevel,celldiameter,anisotropy,numcells,ROIpoints*downsample,ROI,ROIxmin*downsample,ROIymin*downsample);

%create GFP sum and area for each slice
diskSet=calculateSliceGFP(Xorig,diskSet);
  e.diskintensity=diskSet.xymaximavals;
    e.diskGFPsums=diskSet.GFPsums;
      e.diskArea=diskSet.diskArea;
            e.diskMax=diskSet.diskMax; %as 1/31/2019 I dont think this is
           % used and I'm removing it to test
if(savedata)
  
  

    e.numcells=numcells;
    e.celldiameter=celldiameter;
    e.sigma=sigma;
    e.firstroundpoints=diskSet.centeredxymax(nucleiSet.centerindicies,:);
    e.firstroundmaxima=diskSet.xymaximavals(nucleiSet.centerindicies);
    temp=size(diskSet.centeredxymax(nucleiSet.centerindicies,:));
    e.firstroundlength=temp(1);
    e.lowcoverage=length(find(diskSet.xycoverage(nucleiSet.centerindicies)<13));
end
    e.sliceCenters=diskSet.centeredxymax;
    e.xpositions=diskSet.xpositions;
    e.ypositions=diskSet.ypositions;

%iteratively do extra rounds till nothing is left over results are
%concatenated on
nucleiSet=findOverlookedNuclei(X,nucleiSet,diskSet,anisotropy,celldiameter,numcells);
clear('Xorig');
clear ('X');
%conflict resolution
[mergedNucleiSet,e]=resolveConflicts(nucleiSet,diskSet,e,anisotropy,numcells,celldiameter,savedata);



e.finalaveragepoints=single(mergedNucleiSet.averagedpositions);
e.finalpoints=single(mergedNucleiSet.positions);
e.finaldiams=single(mergedNucleiSet.diameters);
e.finalmaximas=mergedNucleiSet.maximas;
e.diams=diskSet.xydetdiameters(nucleiSet.centerindicies); 
e.merged_sliceindicies=mergedNucleiSet.merged_sliceindicies;
   e.allrange=nucleiSet.range;
    e.alllogodds=nucleiSet.logodds;
    e.aspectratio=mergedNucleiSet.aspectratios;
  
if(savedata)
   
    
    %basic calculation saved info
    e.allpoints=diskSet.xymax(nucleiSet.centerindicies,:);
    e.maximavals=diskSet.xymaximavals(nucleiSet.centerindicies);
    e.allcenters=nucleiSet.centers;
   



%matching related saved info
if(~nodata)
    points=mergedNucleiSet.positions;
    mergeinfo2a=e.mergeinfo;
pairmerges=e.pairmerges;

    %read labeled data
[celldata,cellnames]=readnuclei([nucleidir,embryonumber_ed,'\t',num2str(time,'%03d'),'-nuclei']);
%read undedited sn data
%[Ucelldata,Ucellnames]=readnuclei([nucleidir,embryonumber,'_unedited\t',num2str(time,'%03d'),'-nuclei']);

    threshold=.75*celldiameter; %for identity of detected and labeled nuclei

    celllocations=celldata(:,4:6);
%Ucelllocations=Ucelldata(:,4:6);

%stats from initial 1st round matching matching
celllocations(:,1:2)=celllocations(:,1:2).*downsample;
%Ucelllocations(:,1:2)=Ucelllocations(:,1:2).*downsample;
%remember arguments to sub2ind are y,x,z
diameters=celldata(:,7)*downsample;
filteredpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies(1:e.firstroundlength),:);
%compensate for ROI before matching against saved data set
filteredpointsswitched(:,1)=filteredpointsswitched(:,1)+ROIxmin;
filteredpointsswitched(:,2)=filteredpointsswitched(:,2)+ROIymin;
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
filteredpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies,:);
%compensate for ROI before matching against saved data set
filteredpointsswitched(:,1)=filteredpointsswitched(:,1)+ROIxmin;
filteredpointsswitched(:,2)=filteredpointsswitched(:,2)+ROIymin;
[matches2,matchessr2]=compareDetectionWRadius_3(filteredpointsswitched,celllocations,celldata(:,7)*downsample*.5,1.5,anisotropy);

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
filteredpointsswitched=points;
%compensate for ROI before matching against saved data set
filteredpointsswitched(:,1)=filteredpointsswitched(:,1)+ROIxmin;
filteredpointsswitched(:,2)=filteredpointsswitched(:,2)+ROIymin;

[matches2,matchessr2]=compareDetectionWRadius_3(filteredpointsswitched,celllocations,celldata(:,7)*downsample*.5,1.5,anisotropy);
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
