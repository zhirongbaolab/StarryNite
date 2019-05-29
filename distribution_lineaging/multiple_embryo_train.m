lowdim=false; %whether to create full covariance model or independent model for division scoring function (initial creation) false in default usage.

%do all embryo model training for bifurcation model

if(replacemodel)
    
    %these dont exist if original model was run without classifier (or
    %maybe because there are no valid divisions) so putting in replace if
    %kind of a hack
Divdata=[];
NoDivdata=[];
Tripledata=[];


for lin=1:length(lineages)
    Divdata=[Divdata;allbifurcationinfo(lin).Divdata];
NoDivdata=[NoDivdata;allbifurcationinfo(lin).NoDivdata];
Tripledata=[Tripledata;allbifurcationinfo(lin).Tripledata];
end

div_mean=mean(Divdata);


div_std=cov(Divdata);


nodiv_mean=mean(NoDivdata);
nodiv_std=cov(NoDivdata);


if (lowdim)%learn independent feature instead of cov matrix for div scoring
    div_triple_std=zeros(size(div_triple_mean,2));
    stddev=std(Tripledata);
    for ind=1:length(stddev)
        div_triple_std(ind,ind)=stddev(ind);
    end
else
    div_triple_std=cov(Tripledata);
end

div_triple_mean=mean(Tripledata);




    
    trackingparameters.model.div_mean=div_mean;
    trackingparameters.model.div_std=div_std;
    
    trackingparameters.model.div_triple_mean=div_triple_mean;
    trackingparameters.model.div_triple_std=div_triple_std;
    trackingparameters.model.nodiv_mean=nodiv_mean;
    trackingparameters.model.nodiv_std=nodiv_std;
else
    'non replacing model for div '
    
end



FNback=[];
cleanlybad=[];
realdiv=[];


bestFNBackScore=[];
bestFNBackLength=[];
bestFNForwardScoreD1=[];
bestFNForwardScoreD2=[];
    bestFNForwardLengthD1=[];
    bestFNForwardLengthD2=[];
   minsize=[];  
   
       FullyDivLooking=[];
    FNDivLooking=[];
    ncells=[];
    
    FullyFPLooking=[];
    DirtyFPLooking=[];
    DivFPLooking=  [];
    trulyambigious=[];
    alldaughterdata=[];
    allforwarddata=[];
    allbackdata=[];  
    bestFNBackCorrect=[];

for lin=1:length(lineages)
    removed=allbifurcationinfo(lin).removed;
    confidenceData=allbifurcationinfo(lin).confidenceData;
    splitFNMatchScore=allbifurcationinfo(lin).splitFNMatchScore;
    BifurcationMeasures=allbifurcationinfo(lin).BifurcationMeasures;
    ncells=[ncells;allbifurcationinfo(lin).ncells'];
    
    
    
    bestFNBackScore=removed(:,30);
    bestFNBackLength=removed(:,31);
    bestFNForwardScoreD1=removed(:,32);
    bestFNForwardScoreD2=removed(:,33);
    bestFNForwardLengthD1=removed(:,34);
    bestFNForwardLengthD2=removed(:,35);
       bestFNBackCorrect=[bestFNBackCorrect;removed(:,36)];
    minsize=[minsize;min(removed(:,3)',removed(:,4)')'];
    
    localmin=min(removed(:,3)',removed(:,4)')';
    
    %compile all data for bifurcation training
ncells=[];
minbranchforwardconv=[];
minbranchforwardconv_recursive=[];
minsizedaughterbestscore=zeros(size(localmin));
minsizedaughterbestscorelength=zeros(size(localmin));
minbranchforward_solidlength=[];
minbranchforward_timelength=[];
minbranchforwardhops_recursive=[];
minbranchforwardxy=[];
minbranchforwardz=[];
minbranchforwardgapsize=[];
    for i=1:size(removed,1)
        if(localmin(i)==removed(i,3))
            minsizedaughterbestscorelength(i)=bestFNForwardLengthD1(i);
            minsizedaughterbestscore(i)=bestFNForwardScoreD1(i);
            minbranchforwardconv=[minbranchforwardconv;confidenceData.bestforward1conv(i,:)];
            %  minbranchforwardconv_recursive=[minbranchforwardconv_recursive;confidenceData.rforwardd1_conv(i,:)];
            %      minbranchforwardhops_recursive=[minbranchforwardhops_recursive;confidenceData.rforwardd1_consum(i,:)];
            minbranchforward_solidlength=[minbranchforward_solidlength;confidenceData.rforwardd1_solidlength(i)];
            minbranchforward_timelength=[minbranchforward_timelength;confidenceData.rforwardd1_flength(i)];
            minbranchforwardxy=[minbranchforwardxy;splitFNMatchScore.forwardd1xy(i)];
            minbranchforwardz=[minbranchforwardz;splitFNMatchScore.forwardd1z(i)];
            minbranchforwardgapsize=[minbranchforwardgapsize;splitFNMatchScore.forwardd1gapsize(i)];
        else
            
            %       minbranchforwardconv_recursive=[minbranchforwardconv_recursive;confidenceData.rforwardd2_conv(i,:)];
            %     minbranchforwardhops_recursive=[minbranchforwardhops_recursive;confidenceData.rforwardd2_consum(i,:)];
            minbranchforward_solidlength=[minbranchforward_solidlength;confidenceData.rforwardd2_solidlength(i)];
            minbranchforward_timelength=[minbranchforward_timelength;confidenceData.rforwardd2_flength(i)];
            minbranchforwardz=[minbranchforwardz;splitFNMatchScore.forwardd2z(i)];
            minbranchforwardxy=[minbranchforwardxy;splitFNMatchScore.forwardd2xy(i)];
            minsizedaughterbestscorelength(i)=bestFNForwardLengthD2(i);
            minbranchforwardconv=[minbranchforwardconv;confidenceData.bestforward2conv(i,:)];
            minsizedaughterbestscore(i)=bestFNForwardScoreD2(i);
            minbranchforwardgapsize=[minbranchforwardgapsize;splitFNMatchScore.forwardd2gapsize(i)];
            
        end
    end
    'training to not recognize polar bodies'
    cleanlybad=[cleanlybad;(removed(:,11)+removed(:,12)>=localmin)&localmin<=4];%
    FNback=[FNback;removed(:,23)];
     realdiv=[realdiv;logical(removed(:,8))];
  
      %more compact calculation of 6 categories that dont depend on long/short
    FullyDivLooking=[FullyDivLooking;(removed(:,3)>=4&removed(:,4)>=4)&(~(removed(:,24)>0)&~(removed(:,25)>0))];
    FNDivLooking=[FNDivLooking;removed(:,3)>=4&removed(:,4)>=4&(removed(:,24)>0|removed(:,25)>0)];
    
    %local
    small=localmin<4;
    anysmalllacksforwardFN=(removed(:,3)<4&~(removed(:,26)>0))|(removed(:,4)<4&~(removed(:,27)>0));
    %smallhasforwardFN=(removed(:,3)<4&removed(:,26)>0)|(removed(:,4)<4&removed(:,27)>0);
    hasbackwardFN=removed(:,24)>0|removed(:,25)>0;
    
    
    FullyFPLooking=[FullyFPLooking;small&anysmalllacksforwardFN&~hasbackwardFN];
    DirtyFPLooking=[DirtyFPLooking;small&anysmalllacksforwardFN&hasbackwardFN];
    DivFPLooking=  [DivFPLooking;small&~anysmalllacksforwardFN&~hasbackwardFN];
    trulyambigious=[trulyambigious;small&~anysmalllacksforwardFN&hasbackwardFN];
    
    
    alldaughterdata=[alldaughterdata;BifurcationMeasures,confidenceData.bifconvector,localmin,allbifurcationinfo(lin).ncells'];
    allforwarddata=[allforwarddata;minbranchforwardgapsize,minsizedaughterbestscore,minsizedaughterbestscorelength,...
        minbranchforwardconv,minbranchforwardxy,minbranchforwardz,minbranchforward_solidlength,minbranchforward_timelength];
    allbackdata=[allbackdata;splitFNMatchScore.backgapsize,bestFNBackScore,bestFNBackLength,confidenceData.bestbackconv,splitFNMatchScore.backxy,splitFNMatchScore.backz];
    
    
end
   
    
    FullyFPLooking=logical(FullyFPLooking);
    DirtyFPLooking=logical(DirtyFPLooking);
    DivFPLooking= logical(DivFPLooking);
    trulyambigious=logical(trulyambigious);
 

classtags=zeros(size(trulyambigious)); %other=
classtags(logical(realdiv))=1;   %div
%previous dominant correct or other and not
%classtags(~realdiv&((FNback&bestFNBackCorrect)|(~realdiv&~FNback&~cleanlybad&bestFNBackCorrect)))=2; %FN=
classtags(~realdiv&((FNback|bestFNBackCorrect)))=2; %FN=

classtags(cleanlybad&~realdiv&~FNback)=3; %FP


%'stop point for manual training vector'

%{
%random subset selection
%pause here and select subset
%random permutation of data

r=randperm(size(classtags,1));
r=r(1:min(200,size(r,2)));
classtags=classtags(r);
alldaughterdata=alldaughterdata(r,:);
allbackdata=allbackdata(r,:);
allforwarddata=allforwarddata(r,:);
trulyambigious=trulyambigious(r,:);
FullyDivLooking=FullyDivLooking(r,:);
FNDivLooking=FNDivLooking(r,:);
FullyFPLooking=FullyFPLooking(r,:);
DirtyFPLooking=DirtyFPLooking(r,:);
DivFPLooking=DivFPLooking(r,:);
bifurcationPoints=bifurcationPoints(r,:);
bifurcationPoints(:,2)=bifurcationPoints(:,2)+ROIxmin;
bifurcationPoints(:,3)=bifurcationPoints(:,3)+ROIymin;
bifurcationPoints(:,5)=bifurcationPoints(:,5)+ROIxmin;
bifurcationPoints(:,6)=bifurcationPoints(:,6)+ROIymin;
bifurcationPoints(:,8)=bifurcationPoints(:,8)+ROIxmin;
bifurcationPoints(:,9)=bifurcationPoints(:,9)+ROIymin;
bifurcationPoints(:,12)=bifurcationPoints(:,12)+ROIxmin;
bifurcationPoints(:,13)=bifurcationPoints(:,13)+ROIymin


%replace classtags with manual tags and ready to roll

classtags=classtags(1:100);
caseskeep=logical(caseskeep);
classtags=classtags(caseskeep);
alldaughterdata=alldaughterdata(caseskeep,:);
allbackdata=allbackdata(caseskeep,:);
allforwarddata=allforwarddata(caseskeep,:);
trulyambigious=trulyambigious(caseskeep,:);
FullyDivLooking=FullyDivLooking(caseskeep,:);
FNDivLooking=FNDivLooking(caseskeep,:);
FullyFPLooking=FullyFPLooking(caseskeep,:);
DirtyFPLooking=DirtyFPLooking(caseskeep,:);
DivFPLooking=DivFPLooking(caseskeep,:);

trackingparameters.model.div_mean=mean(alldaughterdata(classtags==1,1:2));
trackingparameters.model.div_std=std(alldaughterdata(classtags==1,1:2));

trackingparameters.model.div_triple_mean=mean(alldaughterdata(classtags==1,5:14));
trackingparameters.model.div_triple_std=std(alldaughterdata(classtags==1,5:14));
%}


daughterkeep=logical(ones(size(alldaughterdata,2),1));
backkeep=logical(ones(size(allbackdata,2),1));
forwardkeep=logical(ones(size(allforwarddata,2),1));

for i=1:size(alldaughterdata,2)
    [e1,correctflags,computedlabels,model1]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,daughterkeep)...
        ,allbackdata,allforwarddata(:,forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);
    daughterkeep(i)=0;
    [e2,correctflags,computedlabels,model2]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,daughterkeep)...
        ,allbackdata,allforwarddata(:,forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);
    if(e2>e1)
        daughterkeep(i)=1;
    end
end
for i=1:size(allbackdata,2)
    [e1,correctflags,computedlabels,model1]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,daughterkeep)...
        ,allbackdata(:,backkeep),allforwarddata(:,forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);
    backkeep(i)=0;
    [e2,correctflags,computedlabels,model2]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,daughterkeep)...
        ,allbackdata(:,backkeep),allforwarddata(:,forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);
    if(e2>e1)
    end
end

for i=1:size(allforwarddata,2)
    [e1,correctflags,computedlabels,model1]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,daughterkeep)...
        ,allbackdata(:,backkeep),allforwarddata(:,forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);
    forwardkeep(i)=0;
    [e2,correctflags,computedlabels,model2]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,daughterkeep)...
        ,allbackdata(:,backkeep),allforwarddata(:,forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);
    if(e2>e1)
        forwardkeep(i)=1;
        trackingparameters.bifurcationclassifier=model1;
    else
        trackingparameters.bifurcationclassifier=model2;
    end
end
trackingparameters.bifurcationclassifier.forwardkeep=forwardkeep;
trackingparameters.bifurcationclassifier.backkeep=backkeep;
trackingparameters.bifurcationclassifier.daughterkeep=daughterkeep;

%final

[e1,correctflags,computedlabels,model1]=compute_singlemodel_classificationerror(classtags,alldaughterdata(:,trackingparameters.bifurcationclassifier.daughterkeep)...
        ,allbackdata(:,trackingparameters.bifurcationclassifier.backkeep),allforwarddata(:,trackingparameters.bifurcationclassifier.forwardkeep),trulyambigious,FullyDivLooking,...
        FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking);

