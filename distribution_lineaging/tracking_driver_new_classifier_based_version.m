%lineage driver
%itializes structures
%does easy, gathers candidates, noneasy 1-1, divisions, division classify
%{
if (~isfield(trackingparameters, 'useAveragePoint'))
    trackingparameters.useAveragePoint=false;
end
if(trackingparameters.useAveragePoint)
    for t=1:length(esequence)
        esequence{t}.finalpoints(:,3)=esequence{t}.finalaveragepoints(:,3);
    end
end
%}
%#function distanceCostFunction divScoreModelCostFunction
%#function NaiveBayes
%#function NaiveBayes\predict
%#function predict
%#function wibble
[trackingparameters,esequence]=initializeTrackingStructures(esequence,trackingparameters);

%polar body filter before start tracking


%linking
esequence=linkEasyCases(esequence,trackingparameters);

if (isfield(trackingparameters,'polarbodyfilter')&&trackingparameters.polarbodyfilter)
    windowsize=10;
    maxtime=min(trackingparameters.polarendtime,trackingparameters.endtime);
    for t=1:windowsize:maxtime-windowsize
        sizes=[];
        amaxGFP=[];
        for i=t:t+windowsize
            sizes=[sizes;esequence{i}.finaldiams];
         %   amaxGFP=[amaxGFP;esequence{i}.mdiskMax];
        end
        sizethresh=mean(sizes);
        brightthresh=trackingparameters.PolarThreshold;%mean(amaxGFP)*1.075;
        brightthreshhigh=trackingparameters.PolarThresholdHigh;
    
        if(t>trackingparameters.PolarThreshold2time)
                brightthresh=trackingparameters.PolarThreshold2;%mean(amaxGFP)*1.075;
                brightthreshhigh=trackingparameters.PolarThreshold2High;
        end
        [ esequence] = polarBodyFilter( esequence,t,t+windowsize,sizethresh,brightthresh,brightthreshhigh );
    end
    [ esequence] = polarBodyFilter( esequence,t,maxtime,sizethresh,brightthresh,brightthreshhigh );
    
end
%gather candidates after deleting polar bodies this function is now deleted
%safe
esequence=gatherEndCandidates(esequence,trackingparameters);

%do non gap
%these arent really parameters but internal configuration of endmatch in
%next block
trackingparameters.dotrack=true;
trackingparameters.trackdiv=false;
trackingparameters.tracknondiv=true;
trackingparameters.gapthresh=0;
global answers
answers=[];
wrong=[];
numanswered=[];
result=[];
firstloop=true;
global alldatanondiv;

for offset=trackingparameters.minnondivscore:trackingparameters.nondivscorestep:trackingparameters.maxnondivscore%.875%1 %.1:.0125:.9
    
    nondivthresh=(offset);%.01*1.175^offset;
    alldatadiv=[];
    alldatanondiv=[];
    trackingparameters.endscorethresh_nondiv=nondivthresh;
    esequence=greedyEndScore(esequence,trackingparameters);
   %{
    if(firstloop)
        firstloop=false;
        if(trackingparameters.trainingmode)
            allbifurcationinfo(lin).initialnondivdistribution=alldatanondiv(:,3);
        end
        %trackingparameters.endnodivthreshold=prctile(alldatanondiv(:,3),91);
        %nondivthress(lin)=trackingparameters.maxnondivscore;
    end
    %}
    %compute correct, wrong if in recordanswer mode with answer key
    if(~isempty(answers))
        isdivision=answers(:,15)~=-1;%is chosen answer division
        wrong=[wrong;length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))];
        numanswered=[numanswered;length(find(isdivision)),length(find(~isdivision))];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             [length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))]
        result=[result;[length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))]];
    end
    
end

if(isfield(trackingparameters,'hysteresis')&&trackingparameters.hysteresis)
    esequence=cleanUnlinkedHysteresis(esequence,trackingparameters);  
end

%division closing
trackingparameters.do_enddiv=true;
trackingparameters.dotrack=true;
trackingparameters.trackdiv=true;
trackingparameters.tracknondiv=false;
trackingparameters.gapthresh=0;
answers=[];
%do reasonableish divisions in score order
for offset=trackingparameters.mindivscore:trackingparameters.divscorestep:trackingparameters.maxdivscore %1.5 %.1:.0125:.9
    alldatadiv=[];
    alldatanondiv=[];
    divthresh=offset;%.01*1.175^offset;
      trackingparameters.endscorethresh_div=divthresh;
    esequence=greedyEndScore(esequence,trackingparameters);
    %min score div, min score nondiv,
    %divright,nondivright,answerpresent,answer, answertime, i,t
    %min division is not min
    %isdivision=answers(:,7)~=-1;%is real answer division
    if(~isempty(answers))
        isdivision=answers(:,15)~=-1;%is chosen answer division
        wrong=[wrong;length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))];
        numanswered=[numanswered;length(find(isdivision)),length(find(~isdivision))];
        
        [length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
            ,length(find(~answers(:,4)&~isdivision))]
    end
end



if(~exist('skipbifurcation')||~skipbifurcation)
    
    %do all remaining possible divisions
trackingparameters.endscorethresh_div=inf;
esequence=greedyEndScore(esequence,trackingparameters);
%compute correct/wrong if have answer key
%min score div, min score nondiv,
%divright,nondivright,answerpresent,answer, answertime, i,t
%min division is not min
if(~isempty(answers))
    isdivision=answers(:,15)~=-1;%is chosen answer division
    wrong=[wrong;length(find(~answers(:,3)&isdivision))...
        ,length(find(~answers(:,4)&~isdivision))];
    numanswered=[numanswered;length(find(isdivision)),length(find(~isdivision))];
    [length(find(isdivision)),length(find(~isdivision)),length(find(~answers(:,3)&isdivision))...
        ,length(find(~answers(:,4)&~isdivision))]
end


%test forced error
if(evalforced)
    mkdir(outputdirectory);
    saveGreedyNucleiFiles(esequence,endtime-11,outputdirectory,anisotropy,ROIxmin,ROIymin);
    uneddir=outputdirectory;
    eddir=[basedir,lineages{lin},'_edited\nuclei\'];
    test=evaluate_lineage_error(eddir,uneddir,endtime-11,anisotropy);
    test{37}=test{10}+test{21}+test{32};
    test{38}=test{12}+test{23}+test{34};
    test{39}=sum(test{37}+test{38});
    test{40}=sum(test{37});
    test{41}=sum(test{38});
    forcederrors{lin}=test;
end


%time adjustment back for answer key mode
if(trackingparameters.trainingmode==true)
    trackingparameters.recordanswers=true;
    endtime=endtime-11;
    trackingparameters.endtime=trackingparameters.endtime-11;
else
    trackingparameters.recordanswers=false;%true;
end

%initialize data structures for storing bifurcation results
global bifurcationPoints;
bifurcationPoints=[];
global FNtype;
FNtype=[];
global simpleFNcorrect;
simpleFNcorrect=[];
global removed
global removedi;
global BifurcationMeasures;
global bestCandidateInfo;
global endIndicies;
global bestEndCandidateInfo;
global confidenceData;
global splitFNMatchScore;
global  computedclassificationvector;
global refclassificationvector;
global classround;
classround=[];
refclassificationvector=[];
computedclassificationvector=[];
splitFNMatchScore.forwardd1xy=[];
splitFNMatchScore.forwardd2xy=[];
splitFNMatchScore.forwardd1z=[];
splitFNMatchScore.forwardd2z=[];
splitFNMatchScore.backxy=[];
splitFNMatchScore.backz=[];
splitFNMatchScore.forwardd1gapsize=[];
splitFNMatchScore.forwardd2gapsize=[];
splitFNMatchScore.backgapsize=[];

confidenceData=[];
confidenceData.bifcon=[];
confidenceData.bifconvector=[];
confidenceData.bestbackconv=[];
confidenceData.bestforward1conv=[];
confidenceData.bestforward2conv=[];
confidenceData.rforwardd1_conv=[];
confidenceData.rforwardd1_consum=[];
confidenceData.rforwardd1_flength=[];
confidenceData.rforwardd1_solidlength=[];
confidenceData.rforwardd2_conv=[];
confidenceData.rforwardd2_consum=[];
confidenceData.rforwardd2_flength=[];
confidenceData.rforwardd2_solidlength=[];
endIndicies=[];
bestEndCandidateInfo=[];
bestCandidateInfo=[];
removedi=[];
removed=[];
BifurcationMeasures=[];
esequence=greedydeleteFPbranches(esequence,trackingparameters);

%tally corrected cases
if(trackingparameters.recordanswers)
    %various calculations realated to error correction success for eval
    %purposes.
    minsize=min(removed(:,3)',removed(:,4)');
    cleanlybad=(removed(:,11)+removed(:,12)>=minsize');
    cleanlybad=(cleanlybad&minsize'<6);
    slightlybad=(removed(:,11)+removed(:,12)>=minsize'/2);
    FNback=removed(:,23);
    FNbackwardoption=removed(:,24)>0|removed(:,25)>0;
    realdiv=logical(removed(:,8));
    other=(~slightlybad&~cleanlybad&~removed(:,23)&~realdiv);
    FNforwardfromendofshort=(removed(:,3)<4&removed(:,26)>0)|(removed(:,4)<4&removed(:,27)>0);
    %FNcand exists back but all of them(if more than 1)are small options
    FPbackwardfromFNback=(removed(:,24)>0&removed(:,28)<4)|(removed(:,25)>0&removed(:,29)<4);
    numberofcells=[];
    for i=1:size(removed,1)
        numberofcells=[numberofcells;length(esequence{removed(i,1)}.FP)];
    end
    early=numberofcells<194;
    middle=numberofcells>=194&numberofcells<350;
    rest=~(middle|early);
    all=ones(size(early));
    %global cases
    
    caseofinterest=all;
    
    bestFNBackCorrect=removed(:,36);
    
    
    expected_corrections=[];
    expected_corrections.other=length(find(~FNback&~cleanlybad&~realdiv));
    expected_corrections.FPdet=sum(minsize(cleanlybad&caseofinterest));
    expected_corrections.FPtrack=sum(minsize(cleanlybad&caseofinterest));
    expected_corrections.FNdet=sum(splitFNMatchScore.backgapsize(((simpleFNcorrect==1)&~realdiv&~cleanlybad&caseofinterest))-1);
    %expected_corrections.FPtrack_FN=length(FNback&~cleanlybad&~realdiv);
    expected_corrections.FPtrack_FN=length(find((simpleFNcorrect==1)&~realdiv&~cleanlybad&caseofinterest));
    expected_corrections.FNtrack=sum(splitFNMatchScore.backgapsize(((simpleFNcorrect==1)&~realdiv&~cleanlybad&caseofinterest)));
    expected_corrections.FNdet_swap=sum(splitFNMatchScore.backgapsize((FNtype(:,2)==1&bestFNBackCorrect&~(simpleFNcorrect==1)&~realdiv&~cleanlybad&caseofinterest))-1);
    expected_corrections.FNtrack_swap=sum(splitFNMatchScore.backgapsize((FNtype(:,2)==1&bestFNBackCorrect&~(simpleFNcorrect==1)&~realdiv&~cleanlybad&caseofinterest)));
    expected_corrections.sum=expected_corrections.other+expected_corrections.FPdet+expected_corrections.FPtrack...
        +expected_corrections.FNdet+expected_corrections.FNtrack+expected_corrections.FNdet_swap...
        +expected_corrections.FNtrack_swap+expected_corrections.FPtrack_FN;
end

end


