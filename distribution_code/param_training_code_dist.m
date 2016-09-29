


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

elist=ones(size(esequence));
for example=1:length(esequence)

    ec=esequence{example};
    
    fn=size(ec.finalFN);
    results(tlist(example),elist(example),1)=fn(1);
    fp=size(ec.finalFP);
    results(tlist(example),elist(example),2)=fp(1);
    
    
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





%new ones for current retraining **
starttimes=[1,25,80,181,251,351,501];
endtimes=[24,79,180,250,350,500,2000];
resultstrings={};
counter=1;
for i=1:length(starttimes)

FP1trange=find(FPpoints1nc>=starttimes(i)&FPpoints1nc<endtimes(i));
FN1trange=find(FNpoints1nc>=starttimes(i)&FNpoints1nc<endtimes(i));
FN2trange=find(FNpoints2nc>=starttimes(i)&FNpoints2nc<endtimes(i));
FP2trange=find(FPpoints2nc>=starttimes(i)&FPpoints2nc<endtimes(i));
FP2trangeall=find(FPpoints2allnc>=starttimes(i)&FPpoints2allnc<endtimes(i));
FN3trange=find(FNpoints3nc>=starttimes(i)&FNpoints3nc<endtimes(i));
FP3trange=find(FPpoints3nc>=starttimes(i)&FPpoints3nc<endtimes(i));


['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FN ',num2str(length(FN1trange)),' ',num2str(length(FN2trange)),' ',num2str(length(FN3trange))]
['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FP ',num2str(length(FP1trange)),' ',num2str(length(FP2trangeall)),' ',num2str(length(FP3trange))]
length(find(TPpointsnc>=starttimes(i)&TPpointsnc<endtimes(i)))+length(find(FNpoints2nc>=starttimes(i)&FNpoints2nc<endtimes(i)))

resultstrings{counter}=['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FN ',num2str(length(FN1trange)),' ',num2str(length(FN2trange)),' ',num2str(length(FN3trange))];
counter=counter+1;
resultstrings{counter}=['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FP ',num2str(length(FP1trange)),' ',num2str(length(FP2trangeall)),' ',num2str(length(FP3trange))];
counter=counter+1;
resultstrings{counter}=length(find(TPpointsnc>=starttimes(i)&TPpointsnc<endtimes(i)))+length(find(FNpoints2nc>=starttimes(i)&FNpoints2nc<endtimes(i)));
counter=counter+1;

end
resultstrings=resultstrings';



%}

% pull out all info for parameter tuning
pposition=[];
pintensity=[];
pintensity_e=[];
pnndistance=[];
par=[];
par_e=[];
pcells=[];
pmatching=[];
pstage=[];%1st or 2nd round point

oposition=[];
omerge=[];
osplit=[];
oar=[];
oar_e=[];
ocells=[];
ostage=[];%1st or 2nd round point
omatching=[];
onndist=[];
omademerges=[];%actual live merges made, why are they not same as predicted below
omergear=[];
for example=1:length(elist)
    ec=esequence{example};
  
    pintensity=[pintensity;ec.maximavals];
    if (isfield(ec,'predictedintensity'))
        pintensity_e=[pintensity_e;ec.predictedintensity];
    else
         pintensity_e=[pintensity_e;zeros(length(ec.allcenters),3)];   
    end
    
    pcells=[pcells;ones(length(ec.allcenters),1)*length(ec.allcenters)];
    pmatching=[pmatching;ec.fmatches'];  
    
    stage=[ones(ec.firstroundlength,1);2*ones(length(ec.allcenters)-ec.firstroundlength,1)];
    pstage=[pstage;stage];
   
    if(isfield(ec,'predictedar'))
     par_e=[par_e;ec.predictedar];
    else
         par_e=[par_e;zeros(length(ec.allcenters),3)];
    end
   
    par=[par;ec.allaspectratio];
    points=ec.allpoints;
    spoints=[points(:,1),points(:,2),points(:,3)*anisotropy];
    distancevals=distance(spoints',spoints');
    for j=1:max(size(distancevals))
            distancevals(j,j)=9999;
    end        
   [dist2,I]=min(distancevals);
    dist2=dist2./ec.celldiameter;
      %  dist2=dist2./mean(ec.diams);
    pnndistance=[pnndistance;dist2'];
    pposition=[pposition;ec.allpoints];
%not debugged below    
    %if(isfield(ec,'goodmerges'))
        
   
    %end
    s=size(ec.mergeinfo);
    if(s(1)>0)
        omergear=[omergear;ec.mergear'];
       % omademerges=[omademerges;ec.goodmerges];
    omerge=[omerge;ec.mergeinfo(:,4)];
    osplit=[osplit;ec.mergeinfo(:,3)];
        
    oar=[oar;min([ec.allaspectratio(ec.mergeinfo(:,1))';ec.allaspectratio(ec.mergeinfo(:,2))'])'];
    if(isfield(ec,'predictedar'))
    oar_e=[oar_e;min([ec.predictedar(ec.mergeinfo(:,1))';ec.predictedar(ec.mergeinfo(:,2))'])']; 
    else
        oar_e=[oar_e;zeros(length(ec.mergeinfo),1)];
    end
%    ocells=[ocells;ones(s(1),1)*length(ec.allcenters)];
    if (tlist(example)==firsttimestep) 
        cellnum=firsttimestepnumcells;
    else
       cellnum=length(esequence{example-1}.allcenters);
    end
        %ocells=[ocells;ones(s(1),1)*length(ec.allcenters)];
        ocells=[ocells;ones(s(1),1)*cellnum];
    ostage=[ostage;max([stage(ec.mergeinfo(:,1))';stage(ec.mergeinfo(:,2))'])'];
    omatching=[omatching;ec.fmatches(ec.mergeinfo(:,1))',ec.fmatches(ec.mergeinfo(:,2))'];
    %oposition=[oposition;ec.allpoints(ec.mergeinfo(:,1),:)];
    %oposition=[oposition;ec.allpoints(ec.mergeinfo(:,1),:)];
  oposition=[oposition;max([ec.allpoints(ec.mergeinfo(:,1),1)';ec.allpoints(ec.mergeinfo(:,2),1)'])']; 
  
    for j=1:s(1)
        onndist=[onndist;distancevals(ec.mergeinfo(j,1),ec.mergeinfo(j,2))./ec.celldiameter];
    end
    end
end





%merges version of tuning
%for all times automatically tune thresholds within this range

onn_range=[.55,.6];%[0,1];
omergear_range=[.5,.5];%[.5,2];%[1,2];
osplitmergescale_range=[1,1];%[.1,1.2];%[1,1];%
osplit_range=[0,0]%[1,100];
omerge_range=[-300,0];%[0,20];%[-400,-10];
params=[];

for svalue=1:length(starttimes)
    
    osindicies=(ocells>=starttimes(svalue)&ocells<endtimes(svalue));
    %shouldnt be merged
    tp=omatching(:,1)~=-1&omatching(:,2)~=-1;
    %should be merged
    fp=omatching(:,1)==-1|omatching(:,2)==-1;
    firstround=ostage==1;
    secondround=ostage==2;

    length(find(tp&osindicies));
    length(find(fp&osindicies));
    
    %optimize the flat thresholds on merge score independently to e
    %conservative
    results=[];
    for mergethresh=omerge_range(1):5:omerge_range(2)
        for splitthresh=osplit_range(1):osplit_range(2)
        correct=length(find(fp&osindicies&(osplit<splitthresh&omerge>mergethresh)));
        mistakes=length(find(tp&osindicies&(osplit<splitthresh&omerge>mergethresh)));
        results=[results;correct,mistakes,splitthresh,mergethresh];
        end
    end
    %maximize correctlyremoved-2xfalselyremoved
    [m,mini]=max(results(:,1)-results(:,2)*1);
    mini=find((results(:,1)-results(:,2)*1)>=m, 1, 'last' );%latest param set that is at min error
    

    mergethresh=results(mini,4);
    splitthresh=results(mini,3);
   
    %optimize rest combinatorically
    results=[];
    for nnthresh=onn_range(1):.05:onn_range(2)
        for mergearthresh=omergear_range(1):.1:omergear_range(2)
            for mergesplitscale=osplitmergescale_range(1):.1:osplitmergescale_range(2)
                correct=length(find(fp&osindicies&(onndist<nnthresh|omergear<mergearthresh|omerge>osplit*mergesplitscale|(osplit<splitthresh&omerge>mergethresh))));
                mistakes=length(find(tp&osindicies&(onndist<nnthresh|omergear<mergearthresh|omerge>osplit*mergesplitscale|(osplit<splitthresh&omerge>mergethresh))));
                results=[results;correct,mistakes,nnthresh,mergearthresh,mergesplitscale,splitthresh,mergethresh];
            end
        end
    end
    
    
 %   ['best parameters for stage',num2str(svalue)]
      %maximize correctlyremoved-2xfalselyremoved
    [m,mini]=max(results(:,1)-results(:,2)*1);%most removed compensating for wrong
mini=find((results(:,1)-results(:,2)*1)>=m, 1, 'last' );%latest param set that is at min error
    
    %   'correct merge, wrong merge, nn thresh, mergearthresh, mergesplitscale,splitthresh,mergethresh '
['stage ',num2str(svalue)]
 results(mini,:)

 threshresult=[];
  stagep=(pcells>=starttimes(svalue)&pcells<endtimes(svalue));
   
 for threshsearch=1:2:100
      correctthresh=length(find(pintensity(stagep&pmatching==-1)<threshsearch));
     wrongthresh=length(find(pintensity(stagep&pmatching~=-1)<threshsearch));
        threshresult=[threshresult;threshsearch,correctthresh,wrongthresh];
 end
 [m,mini_thresh]=max(threshresult(:,2)-threshresult(:,3)*1);%most removed compensating for wrong
    

 params=[params; [results(mini,:),threshresult(mini_thresh,1)]];
end
%spit out params in format for pasting into file


nnp='parameters.nndist_merge=[';
mlowp='parameters.mergelower=[';
armergep='parameters.armerge=[';
mergesplitp='parameters.mergesplit=[';
splitp='parameters.split=[';
threshp='parameters.intensitythreshold=[';
for i=1:length(starttimes)
    nnp=[nnp,num2str(params(i,3)),','];
    armergep=[armergep,num2str(params(i,4)),','];
    mergesplitp=[mergesplitp,num2str(params(i,5)),','];
    splitp=[splitp,num2str(params(i,6)),','];
    mlowp=[mlowp,num2str(params(i,7)),','];
    threshp=[threshp,num2str(params(i,8)),','];
 end
    nnp=[nnp,'];']
    mlowp=[mlowp,'];']
    armergep=[armergep,'];']
    mergesplitp=[mergesplitp,'];']
    splitp=[splitp,'];']
    threshp=[threshp,'];']
    
    
return

%just loose 8 by droppping it from zone 4


%optimized for t 6 new images times are not synced btw old and new
length(find(tp&osindicies&(onndist<.6|omergear<1|omerge>osplit|(osplit<12&omerge>.25))))
length(find(fp&osindicies&(onndist<.6|omergear<1|omerge>osplit|(osplit<12&omerge>.25))))

length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<12&omerge>.25))))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<12&omerge>.25))))




%optimized time 5 %for 6 settings
length(find(tp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<0&omerge>-75))&oposition(:,1)>200))
length(find(fp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<0&omerge>-75))&oposition(:,1)>200))
length(find(tp&osindicies&(onndist<.4|omergear<.5|omerge>osplit*1.5|(osplit<0&omerge>-75))&oposition(:,1)<=200))
length(find(fp&osindicies&(onndist<.4|omergear<.5|omerge>osplit*1.5|(osplit<0&omerge>-75))&oposition(:,1)<=200))

%optimized for t 5 new images
length(find(tp&osindicies&(onndist<.6|omergear<1.4|omerge>osplit|(osplit<17&omerge>-25))))
length(find(fp&osindicies&(onndist<.6|omergear<1.4|omerge>osplit|(osplit<17&omerge>-25))))

length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<17&omerge>-25))))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<17&omerge>-25))))



%optimized for t 4 new images
length(find(tp&osindicies&(onndist<.5|omergear<1.3|omerge>osplit|(osplit<20&omerge>-25))))
length(find(fp&osindicies&(onndist<.5|omergear<1.3|omerge>osplit|(osplit<20&omerge>-25))))

length(find(tp&osindicies&(onndist<.5|omerge>osplit|(osplit<20&omerge>-25))))
length(find(fp&osindicies&(onndist<.5|omerge>osplit|(osplit<20&omerge>-25))))



%optimized for t3 8 good 1 fn 
length(find(tp&osindicies&(onndist<.5|omergear<1.2|omerge>osplit*1.5|(osplit<3&omerge>-75))&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.5|omergear<1.2|omerge>osplit*1.5|(osplit<3&omerge>-75))&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.5|omergear<1.2|omerge>osplit*1.5|(osplit<3&omerge>-75))&oposition(:,1)<=225))
length(find(fp&osindicies&(onndist<.5|omergear<1.2|omerge>osplit*1.5|(osplit<3&omerge>-75))&oposition(:,1)<=225))


length(find(tp&osindicies&(omergear<1.2)&oposition(:,1)>225))
length(find(fp&osindicies&(omergear<1.2)&oposition(:,1)>225))
length(find(tp&osindicies&(omergear<1.2)&oposition(:,1)<=225))
length(find(fp&osindicies&(omergear<1.2)&oposition(:,1)<=225))

%optimized for t 3 new images
length(find(tp&osindicies&(onndist<.6|omergear<1.6|omerge>osplit|(osplit<20&omerge>-100))))
length(find(fp&osindicies&(onndist<.6|omergear<1.6|omerge>osplit|(osplit<20&omerge>-100))))

length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<20&omerge>-100))))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<20&omerge>-100))))


%optimized for t2 
%confirm no spatial bias 11good 0fn  17 0fn with .6
%no dif when add merge split
length(find(tp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)<=225))
length(find(fp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)<=225))

length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)<=225))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<8&omerge>-75))&oposition(:,1)<=225))

%optimized for t 2 new images
length(find(tp&osindicies&(onndist<.6|omergear<2|omerge>osplit|(osplit<100&omerge>-200))))
length(find(fp&osindicies&(onndist<.6|omergear<2|omerge>osplit|(osplit<100&omerge>-200))))

length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<100&omerge>-200))))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<100&omerge>-200))))


%optimized for t1 
%7good 0fn %1.5 ar for previous doubled .5 distance
length(find(tp&osindicies&(onndist<.6|omergear<1.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.6|omergear<1.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.6|omergear<1.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)<=225))
length(find(fp&osindicies&(onndist<.6|omergear<1.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)<=225))

length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)<=225))
length(find(fp&osindicies&(onndist<.6|omerge>osplit|(osplit<13&omerge>-75))&oposition(:,1)<=225))

%optimized for t 1 new images
length(find(tp&osindicies&(onndist<.8|omergear<1.6|omerge>osplit|(osplit<100&omerge>-300))))
length(find(fp&osindicies&(onndist<.8|omergear<1.6|omerge>osplit|(osplit<100&omerge>-300))))

length(find(tp&osindicies&(onndist<.8|omerge>osplit|(osplit<100&omerge>-300))))
length(find(fp&osindicies&(onndist<.8|omerge>osplit|(osplit<100&omerge>-300))))





length(find(tp&osindicies&oar<.55&oposition(:,1)>225))
length(find(fp&osindicies&oar<.55&oposition(:,1)>225))
length(find(tp&osindicies&oar<.45))
length(find(fp&osindicies&oar<.45))

length(find(tp&osindicies&oar./oar_e<.75&oposition(:,1)>225))
length(find(fp&osindicies&oar./oar_e<.75&oposition(:,1)>225))
length(find(tp&osindicies&oar./oar_e<.75))
length(find(fp&osindicies&oar./oar_e<.75))


length(find(tp&osindicies&(oar./oar_e<.75|oar<.6)&oposition(:,1)>225))
length(find(fp&osindicies&(oar./oar_e<.75|oar<.6)&oposition(:,1)>225))
length(find(tp&osindicies&(oar./oar_e<.5|oar<.5)))
length(find(fp&osindicies&(oar./oar_e<.5|oar<.5)))


length(find(tp&osindicies&(onndist<.75)&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.75)&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.75)))
length(find(fp&osindicies&(onndist<.75)))
length(find(fp&osindicies))

length(find(tp&osindicies&(onndist<.75|oar<.6)&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.75|oar<.6)&oposition(:,1)>225))
length(find(tp&osindicies&(onndist<.5|oar<.5)))
length(find(fp&osindicies&(onndist<.5|oar<.5)))

length(find(tp&osindicies&(onndist<.5|oar<.5)&oposition(:,1)>225))
length(find(fp&osindicies&(onndist<.5|oar<.5)&oposition(:,1)>225))


return

