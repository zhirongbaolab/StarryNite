


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
for example=1:length(esequence)-1

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

return


ncellsfp=[];
ncellsfn=[];
nfp=[];
nfn=[];
for numcells=2:600
    t=find(FNpoints3nc==numcells);
    if (~isempty(t))
        ncellsfn=[ncellsfn;numcells];
        nfn=[nfn;length(t)];
    end
        t=find(FPpoints3nc==numcells);
    if (~isempty(t))
        ncellsfp=[ncellsfp;numcells];
        nfp=[nfp;length(t)];
    end
end

figure
plot(ncellsfp,nfp,'b');
hold on
plot(ncellsfn,nfn,'r');

%new ones for current retraining **
starttimes=[1,25,80,181,251,351];
endtimes=[24,79,180,250,350,600];

for i=1:length(starttimes)

FP1trange=find(FPpoints1nc>=starttimes(i)&FPpoints1nc<endtimes(i));
FN1trange=find(FNpoints1nc>=starttimes(i)&FNpoints1nc<endtimes(i));
FN2trange=find(FNpoints2nc>=starttimes(i)&FNpoints2nc<endtimes(i));
FP2trange=find(FPpoints2nc>=starttimes(i)&FPpoints2nc<endtimes(i));
FP2trangeall=find(FPpoints2allnc>=starttimes(i)&FPpoints2allnc<endtimes(i));
FN3trange=find(FNpoints3nc>=starttimes(i)&FNpoints3nc<endtimes(i));
FP3trange=find(FPpoints3nc>=starttimes(i)&FPpoints3nc<endtimes(i));


%{

figure
%scatter3(TPpoints(selection,1),TPpoints(selection,2),TPpoints(selection,3)*5,'y')
hold on
scatter3(FNpoints1(FN1trange,1),FNpoints1(FN1trange,2),FNpoints1(FN1trange,3)*anisotropy,'c')
scatter3(FPpoints1(FP1trange,1),FPpoints1(FP1trange,2),FPpoints1(FP1trange,3)*anisotropy,'m')
 % axis([0,325,0,220,0,160])
 axis([0,512,0,512,0,118])
axis('equal')
title(['initial errors t',num2str(starttimes(i)),' to ',num2str(endtimes(i))]);


figure
%scatter3(TPpoints(selection,1),TPpoints(selection,2),TPpoints(selection,3)*anisotropy,'y')
hold on
scatter3(FNpoints2(FN2trange,1),FNpoints2(FN2trange,2),FNpoints2(FN2trange,3)*anisotropy,'c')
%stored is just additional 2nd round ones plot both for moment
 scatter3(FPpoints2(FP2trange,1),FPpoints2(FP2trange,2),FPpoints2(FP2trange,3)*anisotropy,'m')
 %scatter3(FPpoints1(FP1trange,1),FPpoints1(FP1trange,2),FPpoints1(FP1trange,3)*anisotropy,'m')
 %axis([0,325,0,220,0,160])
  axis([0,512,0,512,0,118])
axis('equal')
title(['post 2nd errors t',num2str(starttimes(i)),' to ',num2str(endtimes(i))]);


figure
%scatter3(TPpoints(selection,1),TPpoints(selection,2),TPpoints(selection,3)*anisotropy,'y')
hold on
scatter3(FNpoints3(FN3trange,1),FNpoints3(FN3trange,2),FNpoints3(FN3trange,3)*anisotropy,'c')
scatter3(FPpoints3(FP3trange,1),FPpoints3(FP3trange,2),FPpoints3(FP3trange,3)*anisotropy,'m')
 % axis([0,325,0,220,0,160])
   axis([0,512,0,512,0,118])
axis('equal')
title(['post judgment errors t',num2str(starttimes(i)),' to ',num2str(endtimes(i))]);
%}

%['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FN ',num2str(length(FN1trange)),' ',num2str(length(FN2trange)),' ',num2str(length(FN3trange))]
%['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FP ',num2str(length(FP1trange)),' ',num2str(length(FP2trange)+length(FP1trange)),' ',num2str(length(FP3trange))]
['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FN ',num2str(length(FN1trange)),' ',num2str(length(FN2trange)),' ',num2str(length(FN3trange))]
['time ',num2str(starttimes(i)),' to ',num2str(endtimes(i)),' FP ',num2str(length(FP1trange)),' ',num2str(length(FP2trangeall)),' ',num2str(length(FP3trange))]


length(find(TPpointsnc>=starttimes(i)&TPpointsnc<endtimes(i)))+length(find(FNpoints2nc>=starttimes(i)&FNpoints2nc<endtimes(i)))

end

return



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
    oposition=[oposition;ec.allpoints(ec.mergeinfo(:,1),:)];
 % oposition=[oposition;max([ec.allpoints(ec.mergeinfo(:,1),1)';ec.allpoints(ec.mergeinfo(:,2),1)'])']; 
  
    for j=1:s(1)
        onndist=[onndist;distancevals(ec.mergeinfo(j,1),ec.mergeinfo(j,2))./ec.celldiameter];
    end
    end
end





%merges version of tuning
svalue=6;

osindicies=(ocells>=starttimes(svalue)&ocells<endtimes(svalue));
tp=omatching(:,1)~=-1&omatching(:,2)~=-1;
fp=omatching(:,1)==-1|omatching(:,2)==-1;
firstround=ostage==1;
secondround=ostage==2;

figure
scatter(osplit(tp&osindicies),omerge(tp&osindicies),'b')
hold on
scatter(osplit(fp&osindicies),omerge(fp&osindicies),'r')

figure
scatter(osplit(tp&osindicies&secondround),omerge(tp&osindicies&secondround),'b')
hold on
scatter(osplit(fp&osindicies&secondround),omerge(fp&osindicies&secondround),'r')


length(find(tp&osindicies))
length(find(fp&osindicies))


%test original settings
length(find(tp&osindicies&(onndist<.4|omergear<.55|omerge>1*osplit|(osplit<5&omerge>-20))))
length(find(fp&osindicies&(onndist<.4|omergear<.55|omerge>1*osplit|(osplit<5&omerge>-20))))




%optimized for t 6 new images times are not synced btw old and new
length(find(tp&osindicies&(onndist<.55|omergear<1|omerge>.75*osplit|(osplit<5&omerge>-7))))
length(find(fp&osindicies&(onndist<.55|omergear<1|omerge>.75*osplit|(osplit<5&omerge>-7))))

length(find(tp&osindicies&(onndist<.55|omerge>.75*osplit|(osplit<5&omerge>-20))))
length(find(fp&osindicies&(onndist<.55|omerge>.75*osplit|(osplit<5&omerge>-20))))


.75

.6, 1
alt tuning
length(find(tp&osindicies&(onndist<.59|omergear<1.3|omerge>.15*osplit|(osplit<8&omerge>-50))))
length(find(fp&osindicies&(onndist<.59|omergear<1.3|omerge>.15*osplit|(osplit<8&omerge>-50))))

length(find(tp&osindicies&(onndist<.59|omerge>.15*osplit|(osplit<8&omerge>-50))))
length(find(fp&osindicies&(onndist<.59|omerge>.15*osplit|(osplit<8&omerge>-50))))


for 10 with fudge of 2 edited dist fudge 0
length(find(tp&osindicies&(onndist<.6|omergear<1.1|omerge>1*osplit|(osplit<3.5&omerge>-3))))
length(find(fp&osindicies&(onndist<.6|omergear<1.1|omerge>1*osplit|(osplit<3.5&omerge>-3))))

length(find(tp&osindicies&(onndist<.6|omerge>1.2*osplit|(osplit<3.5&omerge>-3))))
length(find(fp&osindicies&(onndist<.6|omerge>1.2*osplit|(osplit<3.5&omerge>-3))))
edited dist fudge 2

length(find(tp&osindicies&(onndist<.55|omergear<1.1|omerge>.75*osplit|(osplit<3.5&omerge>-5))))
length(find(fp&osindicies&(onndist<.55|omergear<1.1|omerge>.75*osplit|(osplit<3.5&omerge>-5))))

length(find(tp&osindicies&(onndist<.55|omerge>.75*osplit|(osplit<3.5&omerge>-5))))
length(find(fp&osindicies&(onndist<.55|omerge>.75*osplit|(osplit<3.5&omerge>-5))))

edited dist v2 fudge 2 *** best bet for cleaned dist
length(find(tp&osindicies&(onndist<.55|omergear<.5|omerge>.75*osplit|(osplit<3.5&omerge>-10))))
length(find(fp&osindicies&(onndist<.55|omergear<.5|omerge>.75*osplit|(osplit<3.5&omerge>-10))))

length(find(tp&osindicies&(onndist<.55|omerge>.75*osplit|(osplit<3.5&omerge>-10))))
length(find(fp&osindicies&(onndist<.55|omerge>.75*osplit|(osplit<3.5&omerge>-10))))


edited dist v2 fudge 4 *** best bet for cleaned dist
length(find(tp&osindicies&(onndist<.5|omergear<.5|omerge>.75*osplit|(osplit<4.5&omerge>-20))))
length(find(fp&osindicies&(onndist<.5|omergear<.5|omerge>.75*osplit|(osplit<4.5&omerge>-20))))

length(find(tp&osindicies&(onndist<.5|omerge>.75*osplit|(osplit<4.5&omerge>-20))))
length(find(fp&osindicies&(onndist<.5|omerge>.75*osplit|(osplit<4.5&omerge>-20))))



for 10 with fudge of 2 %current best bet
length(find(tp&osindicies&(onndist<.45|omergear<1|omerge>.8*osplit|(osplit<3.5&omerge>-20))))
length(find(fp&osindicies&(onndist<.45|omergear<1|omerge>.8*osplit|(osplit<3.5&omerge>-20))))

length(find(tp&osindicies&(onndist<.45|omerge>.8*osplit|(osplit<3.5&omerge>-20))))
length(find(fp&osindicies&(onndist<.45|omerge>.8*osplit|(osplit<3.5&omerge>-20))))



for dito with alternate 1st round
length(find(tp&osindicies&(onndist<.5|omergear<1|omerge>.5*osplit|(osplit<3.5&omerge>-20))))
length(find(fp&osindicies&(onndist<.5|omergear<1|omerge>.5*osplit|(osplit<3.5&omerge>-20))))

length(find(tp&osindicies&(onndist<.5|omerge>.5*osplit|(osplit<3.5&omerge>-20))))
length(find(fp&osindicies&(onndist<.5|omerge>.5*osplit|(osplit<3.5&omerge>-20))))



for 10 with fudge of 2 and positional bias

length(find(fp&osindicies&(onndist<.45|omergear<1|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)<20))
length(find(fp&osindicies&(onndist<.45|omergear<1|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)>=20))

length(find(tp&osindicies&(onndist<.45|omergear<1|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)<20))
length(find(tp&osindicies&(onndist<.45|omergear<1|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)>=20))



length(find(fp&osindicies&(onndist<.45|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)<20))
length(find(fp&osindicies&(onndist<.55|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)>20))


length(find(tp&osindicies&(onndist<.45|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)<20))
length(find(tp&osindicies&(onndist<.55|omerge>.8*osplit|(osplit<3.5&omerge>-20))&oposition(:,3)>20))



%optimized time 5 %for 6 settings
length(find(tp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<0&omerge>-75))&oposition(:,1)>200))
length(find(fp&osindicies&(onndist<.6|omergear<1.2|omerge>osplit|(osplit<0&omerge>-75))&oposition(:,1)>200))
length(find(tp&osindicies&(onndist<.4|omergear<.5|omerge>osplit*1.5|(osplit<0&omerge>-75))&oposition(:,1)<=200))
length(find(fp&osindicies&(onndist<.4|omergear<.5|omerge>osplit*1.5|(osplit<0&omerge>-75))&oposition(:,1)<=200))

%optimized for t 5 new images
length(find(tp&osindicies&(onndist<.4|omergear<1.1|omerge>.3*osplit|(osplit<17&omerge>-15))))
length(find(fp&osindicies&(onndist<.4|omergear<1.1|omerge>.3*osplit|(osplit<17&omerge>-15))))

length(find(tp&osindicies&(onndist<.4|omerge>.3*osplit|(osplit<17&omerge>-15))))
length(find(fp&osindicies&(onndist<.4|omerge>.3*osplit|(osplit<17&omerge>-15))))
.4 1 no harm
.6 1.4



%optimized for t 4 new images
length(find(tp&osindicies&(onndist<.5|omergear<1.3|omerge>.5*osplit|(osplit<19&omerge>-35))))
length(find(fp&osindicies&(onndist<.5|omergear<1.3|omerge>.5*osplit|(osplit<19&omerge>-35))))

length(find(tp&osindicies&(onndist<.5|omerge>.5*osplit|(osplit<19&omerge>-35))))
length(find(fp&osindicies&(onndist<.5|omerge>.5*osplit|(osplit<19&omerge>-35))))






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
length(find(tp&osindicies&(onndist<.6|omergear<1.6|omerge>.5*osplit|(osplit<20&omerge>-100))))
length(find(fp&osindicies&(onndist<.6|omergear<1.6|omerge>.5*osplit|(osplit<20&omerge>-100))))

length(find(tp&osindicies&(onndist<.6|omerge>.5*osplit|(osplit<20&omerge>-100))))
length(find(fp&osindicies&(onndist<.6|omerge>.5*osplit|(osplit<20&omerge>-100))))


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



%plot mergesplit score components
allgoodmergesv=[];
allbadmergesv=[];
allgoodmergesv2=[];
allbadmergesv2=[];
pall=[];
ebsum=0;
egsum=0;
for example=1:length(elist)
    e=esequence{example};
    p=e.mergeinfo;
   % pall=[p;pall];
    %mratio=p(:,4)./p(:,3);
    mratio=p(:,4);
    matches2=esequence{example}.matches2;
    allbadmergesv=[allbadmergesv;mratio(find(matches2(p(:,1))>0&matches2(p(:,2))>0))];
 allgoodmergesv=[allgoodmergesv;mratio(find(matches2(p(:,1))==-1|matches2(p(:,2))==-1))];
  allbadmergesv2=[allbadmergesv2;p(find(matches2(p(:,1))>0&matches2(p(:,2))>0),3)];
  allgoodmergesv2=[allgoodmergesv2;p(find(matches2(p(:,1))==-1|matches2(p(:,2))==-1),3)];
  
    % should be merged but are odd in that they have negative merge score
    badnet=find((matches2(p(:,1))'==-1|matches2(p(:,2))'==-1)&p(:,4)<-500);
    %sensible merges by our framework
    %badnet=find((matches2(p(:,1))'==-1|matches2(p(:,2))'==-1)&p(:,4)>p(:,3));
    goodnet=find((matches2(p(:,1))'>0&matches2(p(:,2))'>0)&p(:,4)>10);
    goodnet=find((matches2(p(:,1))'>0&matches2(p(:,2))'>0)&p(:,4)<-1000);
    if(length(badnet)>0|length(goodnet)>0)
    ['embryo: ',embryonumbers{elist(example)}, ' ',num2str(elist(example)),' ',num2str(example),' time: ',    num2str(tlist(example))]
    end
      if(length(badnet)>0)
   % 'should be merged but score poorly on merge'
    [e.allpoints(p(badnet,1),:),e.allpoints(p(badnet,2),:),p(badnet,1:4)]
      end
       if(length(goodnet)>0)
   % 'shouldnt be merged but score poorly on split'
  %  [e.allpoints(p(goodnet,1),:),e.allpoints(p(goodnet,2),:),p(goodnet,1:4)]
    end
    ebsum=ebsum+length(badnet);
    egsum=egsum+length(goodnet);
end


figure
scatter(allgoodmergesv2,allgoodmergesv,'r')
hold on
scatter(allbadmergesv2,allbadmergesv,'b')




bm=histc(allbadmergesv,linspace(-200,200,200));
gm=histc(allgoodmergesv,linspace(-200,200,200));
 bar([linspace(-200,200,200);linspace(-200,200,200)]',[bm,gm])
 

%calculate how often windows of logodss do  not overlap completely
 
shortallg=[];
shortallb=[];
nshortallg=0;
nshortallb=0;
shortall=[];
overlaps=[];
r1length=[];
r2length=[];
r1sum=[];
r2sum=[];
r1usum=[];%unique part of r1
r2usum=[];%unique part of r2
r1ovsum=[];%overlapping part of r1
r2ovsum=[];%overlapping part of r2
goodindicies=[];
badindicies=[];
index=1;
for example=1:length(elist)
    e=esequence{example};
    p=e.mergeinfo;
    for i=1:length(p)
        l1=e.alllogodds{p(i,1)};
         l2=e.alllogodds{p(i,2)};
        p1=e.allcenters{p(i,1)};
        p2=e.allcenters{p(i,2)};
        r1=e.allrange{p(i,1)};
        r1length=[r1length,length(r1)];
        r2=e.allrange{p(i,2)};
        r2length=[r2length,length(r2)];
        uniondisks=unique([p1(r1,4);p2(r2,4)]);
        
        %calculate union and unique areas for each 
        r1sum=[r1sum;sum(l1(r1))];
        r2sum=[r2sum;sum(l2(r2))];
        [idisks,IA,IB] = intersect(p1(r1,4),p2(r2,4)); %intersection of disk sets
        l1r=l1(r1);
        l2r=l2(r2);
        r1ovsum=[r1ovsum;sum(l1r(IA))];
        r2ovsum=[r2ovsum;sum(l2r(IB))];
        r1usum=[r1usum,r1sum(index)-r1ovsum(index)];
        r2usum=[r2usum,r2sum(index)-r2ovsum(index)];


        
        overlaps=[overlaps,-1*(length(uniondisks)-length(r1)-length(r2))];

        missing1=0;
        for j=1:length(uniondisks)
            if(length(find(p1(:,4)==uniondisks(j)))==0)
                missing1=missing1+1;
            end
        end
        missing2=0;
        for j=1:length(uniondisks)
            if(length(find(p2(:,4)==uniondisks(j)))==0)
                missing2=missing2+1;
            end
        end

        shortall=[shortall;missing1,missing2];
        matches2=esequence{example}.matches2;
        if((matches2(p(i,1))'==-1|matches2(p(i,2))'==-1))
            goodindicies=[goodindicies,index];
            shortallg=[shortallg;missing1,missing2];
            nshortallg=nshortallg+1;
            
            if(min([r1sum(index),r2sum(index)])>15)
     %       ['embryo: ',embryonumbers{elist(example)}, ' ',num2str(elist(example)),' ',num2str(example),' time: ',    num2str(tlist(example))]
   
    %[e.allpoints(p(i,1),:),e.allpoints(p(i,2),:),p(i,1:4)]
  
        end
            
        end
        if((matches2(p(i,1))>0&matches2(p(i,2))>0))
            
            badindicies=[badindicies,index];
            shortallb=[shortallb;missing1,missing2];
            nshortallb=nshortallb+1;
            
             if(r1sum(index)==r1ovsum(index)|r2sum(index)==r2ovsum(index))
                ['embryo: ',embryonumbers{elist(example)}, ' ',num2str(elist(example)),' ',num2str(example),' time: ',    num2str(tlist(example))]

                [e.allpoints(p(i,1),:),e.allpoints(p(i,2),:),p(i,1:4)]
  
             end
            
        end
        index=index+1;
    end

end

%scatter merge score agains max of nuclei scores
scatter(max([r1sum(goodindicies),r2sum(goodindicies)]'),allgoodmergesv,'r')
hold on
scatter(max([r1sum(badindicies),r2sum(badindicies)]'),allbadmergesv,'b')

scatter(r1ovsum(goodindicies)'./r1sum(goodindicies)',r2ovsum(goodindicies)'./r2sum(goodindicies)','r')
hold on
scatter(r1ovsum(badindicies)'./r1sum(badindicies)',r2ovsum(badindicies)'./r2sum(badindicies)','b')

%scatter(r1ovsum(badindicies)./r1ovsum(badindicies)',r2usum(badindicies)./r2ovsum(badindicies)','b')
[r1ovsum(goodindicies),r1sum(goodindicies),r2ovsum(goodindicies),r2sum(goodindicies)]
[r1ovsum(badindicies),r1sum(badindicies),r2ovsum(badindicies),r2sum(badindicies)]
return



%investigation of temporal properties of size
valid=[];   
merges1=[];%logical indicies into overlap cases
merges2=[];%logical indicies into overlap cases
currentheights=[];   predictedheights=[];
diameters=[];
predicteddiameters=[];  
fragmentmin=[];
fragmentunion=[];
fragmentvalid=[];
fragmentpheight=[];
for time=2:5
    currentpoints=esequence{time}.allpoints;
    previouspoints=esequence{time-1}.allpoints;
    previousheights=[];
 %{   
 %create logical arrays for merge case nuclei
 merges1t=zeros(length(currentpoints),1);
 merges2t=zeros(length(currentpoints),1);
 overlaps=esequence{time}.mergeinfo;
 merges1t(overlaps(:,1))=1;
 merges2t(overlaps(:,2))=1;
 merges1=[merges1,merges1t'];
 merges2=[merges2,merges2t'];
 %}
    %current actual heights
    for i=1:length(currentpoints);
        centerpoints=esequence{time}.allcenters{i};
        centerpoints=centerpoints(esequence{time}.allnewrange{i},3);
        currentheights=[currentheights;max(centerpoints)-min(centerpoints)+1];
      % currentheights=[currentheights;length(esequence{time}.allnewrange{i})];
    end
    %previous actual heights
    for i=1:length(previouspoints);
         centerpoints=esequence{time-1}.allcenters{i};
        centerpoints=centerpoints(esequence{time-1}.allnewrange{i},3);
        previousheights=[previousheights;max(centerpoints)-min(centerpoints)+1];

       %previousheights=[previousheights;length(esequence{time-1}.allnewrange{i})];
    end
    
  %  fragmentmin=[fragmentmin;min(currentheights(overlaps(:,1))*anisotropy./esequence{time}.diams(overlaps(:,1))',currentheights(overlaps(:,2)))*anisotropy./esequence{time}.diams(overlaps(:,2))'];
  %  fragmentvalid=[fragmentvalid,(esequence{time}.fmatches(overlaps(:,1))==-1)|(esequence{time}.fmatches(overlaps(:,2))==-1)];
%{
    
            lpheight=[];
        lpdiams=[];
    % optionally filter previous heights for validity of nucleus
    previousheights=previousheights(find(esequence{time-1}.fmatches~=-1));
    previouspoints=previouspoints(find(esequence{time-1}.fmatches~=-1),:);
    %predict heights
    for i=1:length(currentpoints);
        [minp,maxp,avgp]=predict(currentpoints(i,:),previouspoints,previousheights,4);
        predictedheights=[predictedheights;minp,maxp,avgp];
        [mind,maxd,avgd]=predict(currentpoints(i,:),previouspoints,esequence{time-1}.diams,4);
        predicteddiameters=[predicteddiameters;mind,maxd,avgd];
        lpheight=[lpheight;avgp];
        lpdiams=[lpdiams;avgd];
    end
%}
   valid=logical([valid,logical(esequence{time}.fmatches==-1)]);
    diameters=[diameters,esequence{time}.diams];
 %   fragmentpheight=[fragmentpheight,(lpheight(overlaps(:,1))*anisotropy./lpdiams(overlaps(:,1)))'];
    
end

return

%scatter all points
figure
scatter(predictedheights(valid,3)*anisotropy./predicteddiameters(valid,3),currentheights(valid)*anisotropy./diameters(valid)','r')
hold on
scatter(predictedheights(~valid,3)*anisotropy./predicteddiameters(~valid,3),currentheights(~valid)*anisotropy./diameters(~valid)','b')
axis([0,2,0,2]);

%just nuclei invovled in overlaps
figure
scatter(predictedheights(valid&(merges1|merges2),3)*anisotropy./predicteddiameters(valid&(merges1|merges2),3),currentheights(valid&(merges1|merges2))*anisotropy./diameters(valid&(merges1|merges2))','r')
hold on
scatter(predictedheights(~valid&(merges1|merges2),3)*anisotropy./predicteddiameters(~valid&(merges1|merges2),3),currentheights(~valid&(merges1|merges2))*anisotropy./diameters(~valid&(merges1|merges2))','b')
axis([0,2,0,2]);

%point per merge pair
figure
fragmentvalid=logical(fragmentvalid);
scatter(fragmentpheight(fragmentvalid),fragmentmin(fragmentvalid),'r');
hold on
scatter(fragmentpheight(~fragmentvalid),fragmentmin(~fragmentvalid),'b');

%validpairs=merges1
%minvalues=currentheights(valid&(merges1))*anisotropy./diameters(valid&(merges1))';
%minvalues=min(minvalues,currentheights(valid&(merges2))*anisotropy./diameters(valid&(merges2))');
%minvaluesi=currentheights(~valid&(merges1))*anisotropy./diameters(~valid&(merges1))';
%minvaluesi=min(minvaluesi,currentheights(~valid&(merges2))*anisotropy./diameters(~valid&(merges2))');

%scatter(predictedheights(valid&(merges1),3)*anisotropy./predicteddiameters(valid&(merges1),3),minvalues,'r')
%hold on
%scatter(predictedheights(~valid&(merges1),3)*anisotropy./predicteddiameters(~valid&(merges1),3),minvaluesi','b')
%axis([0,2,0,2]);



scatter(predictedheights(valid,2),currentheights(valid),'r')
hold on
scatter(predictedheights(~valid,2),currentheights(~valid),'b')

return


%num each case w missing overalp
 
 foverlap=overlaps./(r1length+r2length)
 minsize=min([r1length;r2length]);
 foverlaps2=overlaps./min([r1length;r2length]);
 
 
 
 
for example=1:length(elist)
    esequence{example}
end




logoddsgood=[];
logoddsbad=[];

for example=1:length(elist)
    ec=esequence{example};
    
    goodnuc=find(ec.fmatches>0);
    badnuc=find(ec.fmatches==-1);
    logodds=[];
    
    for t=1:length(ec.allcenters)
        %when you implement (for the 3d time) aspect ratio
        %remember logodds less than threshold are resequencey zero
        lo=ec.alllogodds{t};
        losum=sum(lo(ec.allnewrange{t}));
        if(length(ec.allnewrange{t})~=0)
        losum=losum./(length(ec.allnewrange{t})-1);
        %else its zero already
        end
        %{
        if(losum<0)
            t
            example
        end
        %}
        logodds=[logodds,losum];
    end
    
    logoddsgood=[logoddsgood,logodds(goodnuc)];
    logoddsbad=[logoddsbad,logodds(badnuc)];
    
    
   
  
end

maxloval=max([logoddsgood,logoddsbad]);
minloval=min([logoddsgood,logoddsbad]);
minloval=-5
good=histc(logoddsgood,linspace(minloval,maxloval,100));
bad=histc(logoddsbad,linspace(minloval,maxloval,100));
bar([linspace(minloval,maxloval,100);linspace(minloval,maxloval,100)]',[good;bad]');
