clear
load('G:\My Documents\MATLAB\clean_distributions.mat');



%{
allarm=[];
allarmin=[];
allarp=[];
%}
allvalid=[];


suffix='';

elist=[1,2,3,4];
tlist=linspace(2,195,(195-2+1));
firsttimestep=2;
firsttimestepdiam=80;
firsttimestepnumcells=4; %guess at cell stage for parameter set choice




elist=[1,2];
tlist=[190,191,192];
firsttimestep=190;
firsttimestepdiam=17*2;%not used when nodata=false
firsttimestepnumcells=300; %guess at cell stage for parameter set choice

elist=[1,2,3,4];
tlist=linspace(180,195,(195-180+1));
firsttimestep=180;
firsttimestepdiam=34;
firsttimestepnumcells=300; %guess at cell stage for parameter set choice

elist=1;
tlist=185;
firsttimestep=185;
firsttimestepdiam=17*2;%not used when nodata=false
firsttimestepnumcells=300; %guess at cell stage for parameter set choice

newscope=true;
LSM=true;%mouse data
if (newscope)
    nodata=true;%whether to use SN data to match online
    nodatause=true;%whether to use SN data for diameter and stored bottom
    savedata=true;
    singlevolume=false;
    elist=[1];    
    firsttimestep=1;
    if(LSM)
        tlist=linspace(1,20,20); 
        firsttimestepdiam=9;
       
        firsttimestepnumcells=500;
    else

    tlist=linspace(1,450,(450-1+1));
    % tlist=linspace(1,20,(20-1+1));
     firsttimestepnumcells=4;
    firsttimestepdiam=30;
     %guess at cell stage for parameter set choice
    
  %      tlist=200;
    % tlist=linspace(1,20,(20-1+1));
   % firsttimestep=200;
   % firsttimestepdiam=15;
   % firsttimestepnumcells=200;
    end
else
    nodata=false;%whether to use SN data to match online
    nodatause=true;%whether to use SN data for diameter and stored bottom
    savedata=true;
    singlevolume=false;
end

global parameters;

parameters.sigma=1;

if(newscope)
    parameters.staging=[25,102,181,251,351];
   % parameters.intensitythreshold=10;
   if (LSM)
      parameters.intensitythreshold=[5,5,5,5,5,5];
   else
   parameters.intensitythreshold=[10,20,30,30,30,30];
   end
    %parameters.intensitythreshold=.00025;
    parameters.rangethreshold=[100,100,15,8,6,3];
    parameters.nndist_merge=[.6,.6,.6,.5,.5,.4];
    parameters.mergelower=[-75,-75,-75,-75,-30,-75];
    parameters.armerge=[1.6,1.6,1.2,1.2,.5,.5];
    parameters.mergesplit=[1,1,1,1.2,1,1];
    parameters.split=[13,13,8,3,5,0];

else
    parameters.staging=[102,181,251,351];
    parameters.intensitythreshold=.02;
    parameters.rangethreshold=[100,15,8,6,3];
    parameters.nndist_merge=[.6,.6,.5,.5,.4];
    parameters.mergelower=[-75,-75,-75,-30,-75];
    parameters.armerge=[1.6,1.2,1.2,.5,.5];
    parameters.mergesplit=[1,1,1.2,1,1];
    parameters.split=[13,8,3,5,0];

end


SNoutput=true;
%config stuff for files
if (newscope)
    rednuclei=false;
    outputSlice=true;
    downsample=1;
    if(LSM) 
        anisotropy=1.5/1.24*downsample;

        eachslice=[98];
        embryodir='L:\santella\mouse\';


        nucleidir='L:\santella\mouse\'
        embryonumbers={'on.lsm'};
        embryonumbers_ed={'on.lsm'};

    else
        anisotropy=1/.254*downsample; %this is wrong
        eachslice=[30];
        embryodir='L:\newembryosuntiltwitching\originalmetamorphdata\s2\';
        embryodir='L:\duz\project\Imaging\BV24_WT_20090921\';
        nucleidir='L:\newembryosuntiltwitching\originalmetamorphdata\s2\';
        nucleidir='L:\duz\project\Imaging\BV24_WT_20090921\';


        embryonumbers={'BV24_wt_061710_s2'};
        embryonumbers={'ZD_BV24_WT_20090921_s1'};
        embryonumbers_ed={'BV24_wt_061710_s2'};
        embryonumbers_ed={'ZD_BV24_WT_20090921_s1'};
    end
  
  zeropadding=false;
   ROI=false;
ROIxmin=137;
ROIxmax=346;
ROIymin=270;
ROIymax=422; 
   
else
    downsample=.5;   
    anisotropy=10*downsample;
    eachslice=[31;35;31;31];
    embryonumbers = {'082005','081505','081905','090205'};
    embryonumbers_ed = {'082005.2','081505.2','081905.2','090205.2'};
    embryolevels=[25,28,26,29]; %last plane maxima is legal on
    rootembryodir='C:\images\';
    nucleidir='L:\bao\dev\data analysis\data\';

end
%outputmatlabfile=['G:\My Documents\MATLAB\distribution_code\',embryonumbers{1},'_newpipeline_looserthresh.mat'];
parameters.regions=cell(5,1); 
if (~newscope)
%regional definitions of subareas to process differently
parameters.regions=cell(6,1); 
parameters.regions{4}.rangethreshold=16;
parameters.regions{4}.area=[450*downsample,720*downsample,0,512*downsample,0,35];
parameters.regions{4}.armerge=.65;
parameters.regions{4}.split=11;
parameters.regions{4}.nndist_merge=.65;

parameters.regions{5}.rangethreshold=8;
parameters.regions{5}.nndist_merge=.6;
parameters.regions{5}.armerge=1.2;
parameters.regions{5}.area=[450*downsample,720*downsample,0,512*downsample,0,35];
else
parameters.regions=cell(6,1); 
end



eall={};
tic
bottomdata={};
for example=1:length(elist)
    emind=elist(example);
    if(singlevolume&~nodata)
        zlevel=embryolevels(emind);
    else
        zlevel=eachslice(emind);
    end
    
    slices=eachslice(emind);
    embryonumber=embryonumbers{emind};
    nucleibase=[nucleidir,embryonumber,'\'];
    embryonumber_ed=embryonumbers_ed{emind};
    if (~newscope)
    embryodir=[rootembryodir,embryonumber,'\image\tif\'];
    end
    esequence={};
    processSequence;
    eall={eall{:},esequence{:}};
end


toc
%nothing to tally if didnt match online against saved data
if(~nodata&savedata)




    eacc.initialFN=0;
    eacc.initialFP=0;
    eacc.introducedmergeFN= 0;
    eacc.correctmerges=0;
    eacc.secondroundFP= 0;
    eacc.correct2ndroundnuclei=0;
    eacc.zerodiskcorrect=  0;
    eacc.zerodiskwrong= 0;
    eacc.finalFN= 0;
    eacc.finalFP=  0;
    eacc.nuclei=0;

    eaccv.initialFN=[];
    eaccv.initialFP=[];
    eaccv.introducedmergeFN= [];
    eaccv.correctmerges=[];
    eaccv.secondroundFP= [];
    eaccv.correct2ndroundnuclei=[];
    eaccv.zerodiskcorrect= [];
    eaccv.zerodiskwrong= [];
    eaccv.finalFN=[];
    eaccv.finalFP=[];
    eaccv.nuclei=[];
    
        eaccv.secfpsize=[];
        eaccv.sectpsize=[];

        TPpoints=[];
        FPpoints1=[];
        FNpoints1=[];       
        FPpoints2=[];
        FPpoints2all=[];
        FNpoints2=[]; 
        FPpoints3=[];
        FNpoints3=[];
        TPpointsnc=[];
        FPpoints1nc=[];
        FNpoints1nc=[];       
        FPpoints2nc=[];
        FNpoints2nc=[]; 
                FPpoints2allnc=[];
        FPpoints3nc=[];
        FNpoints3nc=[];
        results=zeros(195,4,2);%store all final FN FP
%        for emb=1:4
%            for timepoint=175:195
% example=(timepoint-1)*emb;               
for example=1:length(eall)

    ec=eall{example};
    
    fn=size(ec.finalFN);
   % results(tlist(example),elist(example),1)=fn(1);
    fp=size(ec.finalFP);
   % results(tlist(example),elist(example),2)=fp(1);
    
    
    eacc.initialFN=eacc.initialFN+numel(ec.initialFN)/3;
    eacc.initialFP=eacc.initialFP+numel(ec.initialFP)/3;
    eacc.introducedmergeFN= eacc.introducedmergeFN+numel(ec.introducedmergeFN)/4;
    eacc.correctmerges=eacc.correctmerges+numel(ec.correctmerges)/4;
    eacc.secondroundFP= eacc.secondroundFP+numel(ec.secondroundFP)/3;
    eacc.correct2ndroundnuclei= eacc.correct2ndroundnuclei+numel(ec.correct2ndroundnuclei)/3;
    eacc.zerodiskcorrect=  eacc.zerodiskcorrect+numel(ec.zerodiskcorrect)/3;
    eacc.zerodiskwrong= eacc.zerodiskwrong+numel(ec.zerodiskwrong)/3;
    eacc.finalFN= eacc.finalFN+numel(ec.finalFN)/3;
    eacc.finalFP=  eacc.finalFP+numel(ec.finalFP)/3;
    eacc.nuclei=eacc.nuclei+ec.nuclei;
    
    range=ec.allnewrange;
    for t=1:length(ec.secondroundFPind)
        eaccv.secfpsize=[eaccv.secfpsize;length(range{ec.secondroundFPind(t)})];
    end
    for t=1:length(ec.correct2ndroundnucleiind)
        eaccv.sectpsize=[ eaccv.sectpsize;length(range{ec.correct2ndroundnucleiind(t)})];
    end
    eaccv.initialFN=[eaccv.initialFN;ec.initialFN];
    eaccv.initialFP=[eaccv.initialFP;ec.initialFP];
    eaccv.introducedmergeFN= [eaccv.introducedmergeFN;ec.introducedmergeFN];
    eaccv.correctmerges=[eaccv.correctmerges;ec.correctmerges];
    eaccv.secondroundFP= [eaccv.secondroundFP;ec.secondroundFP];
    eaccv.correct2ndroundnuclei= [eaccv.correct2ndroundnuclei;ec.correct2ndroundnuclei];
    eaccv.zerodiskcorrect=  [eaccv.zerodiskcorrect;ec.zerodiskcorrect];
    eaccv.zerodiskwrong= [eaccv.zerodiskwrong;ec.zerodiskwrong];
    eaccv.finalFN= [eaccv.finalFN;ec.finalFN];
    eaccv.finalFP=  [eaccv.finalFP;ec.finalFP];
    
        TPpoints=[TPpoints;ec.allpoints(find(ec.fmatches~=-1),:)];
        FPpoints1=[FPpoints1;ec.initialFP];
        FNpoints1=[FNpoints1;ec.initialFN];       
        FPpoints2=[FPpoints2;ec.secondroundFP];%just 2nd round ones
        
        FPpoints2all=[FPpoints2all;ec.allpoints(ec.fmatches==-1,:)];
        FNpoints2=[FNpoints2;ec.all2ndroundFN];        
        FPpoints3=[FPpoints3;ec.finalFP];
        FNpoints3=[FNpoints3;ec.finalFN];
        
        TPpointsnc=[TPpointsnc;ones(size(ec.allpoints(find(ec.fmatches~=-1),1)))*(ec.nuclei)];
        FPpoints1nc=[FPpoints1nc;ones(size(ec.initialFP(:,1)))*(ec.nuclei)];
        FNpoints1nc=[FNpoints1nc;ones(size(ec.initialFN(:,1)))*(ec.nuclei)];       
        FPpoints2nc=[FPpoints2nc;ones(size(ec.secondroundFP(:,1)))*(ec.nuclei)];
        FPpoints2allnc=[FPpoints2allnc;ones(length(find(ec.fmatches==-1)),1)*ec.nuclei];
        FNpoints2nc=[FNpoints2nc;ones(size(ec.all2ndroundFN(:,1)))*(ec.nuclei)]; 
        FPpoints3nc=[FPpoints3nc;ones(size(ec.finalFP(:,1)))*(ec.nuclei)]; 
        FNpoints3nc=[FNpoints3nc;ones(size(ec.finalFN(:,1)))*(ec.nuclei)]; 
  
  %          end
        end
eacc



end
%clear X;
%clear Xr;
%save(outputmatlabfile);



