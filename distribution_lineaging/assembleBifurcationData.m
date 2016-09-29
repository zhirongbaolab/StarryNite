function  [d1cand,d2cand,d1candt,d2candt,d1length,d2length,...
    bestFNBackCorrect,bestmatchings,FNbackcand1lengths,FNbackcand2lengths,...
    bestFNForwardLengthD1,bestFNForwardLengthD2,....
    bestmatchingsplayerstart,bestmatchingplayerend,bestdaughter,bestIndex,...
    splitscores_div,alldaughterdata,allforwarddata,allbackdata,count] ....
    = assembleBifurcationData( esequence,t,i,trackingparameters,count )
%traverse and examine bifurcation, add data about it to global variables
%that store it  for training as well as passing back
%variables necessary for classification processing after classificaton

global bifurcationPoints;
global simpleFNcorrect;
global BifurcationMeasures;
global candidateInfo;
global bestCandidateInfo;
global endIndicies;
global bestEndCandidateInfo;
global confidenceData;
global splitFNMatchScore;
global removed;%store answer for post tunign


%normalization value avg nn distance
normv=mean(esequence{t}.selfdistance);
if(isinf(normv))
    normv=mean(esequence{t}.finaldiams);
end

%calculate lengths of all fragments involved in this
d1length=traverse_forward(esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1));
d2length=traverse_forward(esequence,esequence{t}.suc_time(i,2),esequence{t}.suc(i,2));

minlength=min(d1length,d2length);
%compute confidence features for daughters
[d1convector,d1consum]=computeConfidence(esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),minlength);
[d2convector,d2consum]=computeConfidence(esequence, esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),minlength);
confidenceData.bifcon=[confidenceData.bifcon;(d1consum+d2consum)/(minlength*2)];
confidenceData.bifconvector=[confidenceData.bifconvector;mean([d1convector;d2convector])];

d1cur=esequence{t}.suc(i,1);
d2cur=esequence{t}.suc(i,2);
tcur1=esequence{t}.suc_time(i,1);
tcur2=esequence{t}.suc_time(i,2);


localbpoint=[t,esequence{t}.finalpoints(i,:),...
    esequence{tcur1}.finalpoints(d1cur,:),esequence{tcur2}.finalpoints(d2cur,:)];




if(trackingparameters.recordanswers)
    dFNf=esequence{tcur1}.correct_suc_time(d1cur,1)~=-1&esequence{tcur1}.correct_suc_time(d1cur)~=tcur1+1;
    dFNb=esequence{tcur1}.correct_pred_time(d1cur)~=-1&esequence{tcur1}.correct_pred_time(d1cur)~=tcur1-1;
    
    %else
    dFNf=dFNf|esequence{tcur2}.correct_suc_time(d2cur,1)~=-1&esequence{tcur2}.correct_suc_time(d2cur)~=tcur2+1;
    dFNb=dFNb|esequence{tcur2}.correct_pred_time(d2cur)~=-1&esequence{tcur2}.correct_pred_time(d2cur)~=tcur2-1;
else
    dFNf=-1;
    dFNb=-1;
end

%assemble backward candidates
[d1cand,d1candt]=findBackwardCandidatesTime(esequence,d1cur,tcur1,trackingparameters);
[d2cand,d2candt]=findBackwardCandidatesTime(esequence,d2cur,tcur2,trackingparameters);

candidateInfo(count).t=t;
candidateInfo(count).i=i;
candidateInfo(count).d1cand=d1cand;
candidateInfo(count).d2cand=d2cand;
candidateInfo(count).d1candt=d1candt;
candidateInfo(count).d2candt=d2candt;


FNbackcand1=length(d1cand)./trackingparameters.interval;
FNbackcand2=length(d2cand)./trackingparameters.interval;
FNbackcand1lengths=[];
FNbackcand2lengths=[];
for j=1:length(d1cand)
    FNbackcand1lengths=[FNbackcand1lengths;traverse_backdivstop(esequence,d1candt(j),d1cand(j))./trackingparameters.interval];
end
for j=1:length(d2cand)
    FNbackcand2lengths=[FNbackcand2lengths;traverse_backdivstop(esequence,d2candt(j),d2cand(j))./trackingparameters.interval];
end

%compute best FN option backward
[bestdaughter,bestIndex,bestFNBackScore,iscorrect,bestmatchings,bestmatchingsplayerstart,bestmatchingplayerend]=computeBestFNBackOption(d1cand,d1candt,d2cand,d2candt,tcur1,d1cur,d2cur,esequence,trackingparameters);
bestCandidateInfo(count).t=t;
bestCandidateInfo(count).i=i;
if (bestdaughter==-1)
    bestCandidateInfo(count).cand=0;
    bestCandidateInfo(count).candt=0;
    
    confidenceData.bestbackconv=[confidenceData.bestbackconv;zeros(1,6)];
    splitFNMatchScore.backxy=[splitFNMatchScore.backxy;-1];
    splitFNMatchScore.backz=[splitFNMatchScore.backz;-1];
    splitFNMatchScore.backgapsize=[splitFNMatchScore.backgapsize;-1];
    
else
    %if exists compute confidence vector for best back gap target
    if(bestdaughter==1)
        bestCandidateInfo(count).cand=d1cand(bestIndex);
        bestCandidateInfo(count).candt=d1candt(bestIndex);
        
        [bbconvector,bbconsum]=computeConfidenceback(esequence,d1candt(bestIndex),d1cand(bestIndex),FNbackcand1lengths(bestIndex));
        v=bbconvector;
        if(size(v,1)>1)
            v=mean(bbconvector);
        end
        confidenceData.bestbackconv=[confidenceData.bestbackconv;v];
        splitFNMatchScore.backgapsize=[splitFNMatchScore.backgapsize;(tcur1-d1candt(bestIndex))./trackingparameters.interval];
        
        
    else
        bestCandidateInfo(count).cand=d2cand(bestIndex);
        bestCandidateInfo(count).candt=d2candt(bestIndex);
        
        [bbconvector,bbconsum]=computeConfidenceback(esequence,d2candt(bestIndex),d2cand(bestIndex),FNbackcand2lengths(bestIndex));
        v=bbconvector;
        if(size(v,1)>1)
            v=mean(bbconvector);
        end
        confidenceData.bestbackconv=[confidenceData.bestbackconv;v];
        splitFNMatchScore.backgapsize=[splitFNMatchScore.backgapsize;(tcur2-d2candt(bestIndex))./trackingparameters.interval];
        
        
    end
end
if(bestCandidateInfo(count).candt~=0)
bifurcationPoints=[bifurcationPoints;localbpoint,bestCandidateInfo(count).candt,...
    esequence{bestCandidateInfo(count).candt}.finalpoints(bestCandidateInfo(count).cand,:)];
else
    bifurcationPoints=[bifurcationPoints;localbpoint,0,0,0,0];

end

%length of best fn back
bestFNBackCorrect=iscorrect;
if(bestdaughter==-1)
    bestFNBackLength=-1;
else
    if (bestdaughter==1);
        bestFNBackLength=FNbackcand1lengths(bestIndex);
    else
        bestFNBackLength=FNbackcand2lengths(bestIndex);
    end
end

%compute x,y z distances for best match
bestbackind=-1;
if(bestdaughter==1)
    bestdindex=d1cur;
    bestbackind=d1cand(bestIndex);
    bestbacktime=d1candt(bestIndex);
    splitFNMatchScore.backxy=[splitFNMatchScore.backxy;distance(esequence{d1candt(bestIndex)}.finalpoints(d1cand(bestIndex),1:2)',esequence{tcur1}.finalpoints(d1cur,1:2)')...
        /normv];
    splitFNMatchScore.backz=[splitFNMatchScore.backz;abs(esequence{d1candt(bestIndex)}.finalpoints(d1cand(bestIndex),3)-esequence{tcur1}.finalpoints(d1cur,3))...
        *trackingparameters.anisotropyvector(3)/normv];
else
    if(bestdaughter==2)
        bestdindex=d2cur;
        bestbackind=d2cand(bestIndex);
        bestbacktime=d2candt(bestIndex);
        splitFNMatchScore.backxy=[splitFNMatchScore.backxy;distance(esequence{d2candt(bestIndex)}.finalpoints(d2cand(bestIndex),1:2)',esequence{tcur2}.finalpoints(d2cur,1:2)')...
            /normv];
        splitFNMatchScore.backz=[splitFNMatchScore.backz;abs(esequence{d2candt(bestIndex)}.finalpoints(d2cand(bestIndex),3)-esequence{tcur2}.finalpoints(d2cur,3))...
            *trackingparameters.anisotropyvector(3)/normv];
    end
end

%walk forward along 1 till at end and find candidates
for j=1:d1length-1
    d1curbk=d1cur;
    d1cur=esequence{tcur1}.suc(d1cur,1);
    tcur1=esequence{tcur1}.suc_time(d1curbk,1);
end
endIndicies(count).d1end=d1cur;
endIndicies(count).d1endt=tcur1;

[d1endcand,d1endcandt] =findForwardCandidatesTime(esequence,d1cur,tcur1,trackingparameters);

branchendFNoption1=length(d1endcand);
for j=1:d2length-1
    d2curbk=d2cur;
    d2cur=esequence{tcur2}.suc(d2cur,1);
    tcur2=esequence{tcur2}.suc_time(d2curbk,1);
end
endIndicies(count).d2end=d2cur;
endIndicies(count).d2endt=tcur2;

[d2endcand,d2endcandt]=findForwardCandidatesTime(esequence,d2cur,tcur2,trackingparameters);
branchendFNoption2=length(d2endcand);                %{

%compute lengths of forward options of both daughters
d1forwardlengths=[];
d2forwardlengths=[];
for j=1:length(d1endcand);
    d1forwardlengths=[d1forwardlengths,traverse_forward(esequence,d1endcandt(j),d1endcand(j))./trackingparameters.interval];
end
for j=1:length(d2endcand);
    d2forwardlengths=[d2forwardlengths,traverse_forward(esequence,d2endcandt(j),d2endcand(j))./trackingparameters.interval];
end

%compute scores of candidates at ends for both daughters
%(useful if both are short)
if(isempty(d1endcand))
    bestFNForwardScoreD1=-1;
    bestFNForwardLengthD1=-1;
    confidenceData.bestforward1conv=[confidenceData.bestforward1conv;zeros(1,6)];
    confidenceData.rforwardd1_conv=[confidenceData.rforwardd1_conv;zeros(1,6)];
    confidenceData.rforwardd1_consum=[confidenceData.rforwardd1_consum;-1];
    confidenceData.rforwardd1_flength=[confidenceData.rforwardd1_flength;-1];
    confidenceData.rforwardd1_solidlength=[confidenceData.rforwardd1_solidlength;-1];
    splitFNMatchScore.forwardd1xy=[splitFNMatchScore.forwardd1xy;-1];
    splitFNMatchScore.forwardd1z=[splitFNMatchScore.forwardd1z;-1];
    splitFNMatchScore.forwardd1gapsize=[splitFNMatchScore.forwardd1gapsize;-1];
    
else
    [scores_gapd1,matchingsd1,playersstartd1,playersendd1]=gapScore(esequence,d1cur,tcur1,d1endcand,d1endcandt,trackingparameters);
    [bestFNForwardScoreD1,minid1]=min(scores_gapd1);
    bestFNForwardLengthD1=d1forwardlengths(minid1);
    
    splitFNMatchScore.forwardd1gapsize=[splitFNMatchScore.forwardd1gapsize;d1endcandt(minid1)-tcur1];
    splitFNMatchScore.forwardd1xy=[splitFNMatchScore.forwardd1xy;distance(esequence{d1endcandt(minid1)}.finalpoints(d1endcand(minid1),1:2)',esequence{tcur1}.finalpoints(d1cur,1:2)')...
        /normv];
    splitFNMatchScore.forwardd1z=[splitFNMatchScore.forwardd1z;abs(esequence{d1endcandt(minid1)}.finalpoints(d1endcand(minid1),3)-esequence{tcur1}.finalpoints(d1cur,3))...
        *trackingparameters.anisotropyvector(3)/normv];
    
    
    origparent=esequence{d1endcandt(minid1)}.pred(d1endcand(minid1));
    if(origparent~=-1)%option originating in a division
        parentd1=esequence{d1endcandt(minid1)-1}.suc(origparent,1);
        parentd2=esequence{d1endcandt(minid1)-1}.suc(origparent,2);
        %note here lengths are in frames for procesing not
        %minutes as cue
        d1lengthforward=traverse_forward(esequence,d1endcandt(minid1),parentd1);
        d2lengthforward=traverse_forward(esequence,d1endcandt(minid1),parentd2);
        targetbiflength=min(d1lengthforward,d2lengthforward); %overwrite length of target with min length of bifurcation
        bestFNForwardLengthD1=targetbiflength;
        [bbconvector,bbconsum]=computeConfidence(esequence,d1endcandt(minid1),parentd1,bestFNForwardLengthD1);
        v=bbconvector;
        [bbconvector,bbconsum]=computeConfidence(esequence,d1endcandt(minid1),parentd2,bestFNForwardLengthD1);
        
        v=[v;bbconvector];
        if(size(v,1)>1)
            v=mean(v);
        end
        %compute recursive solidity of smaller branch of best
        %candidate pair
        [convectord1,consumd1,flengthd1, solidlengthd1]=computeConfidenceRecursive(esequence,d1endcandt(minid1),d1endcand(minid1),trackingparameters);
        
    else %option not originating in division
        
        %if its not a bifurcation just compute the dangling bit
        [bbconvector,bbconsum]=computeConfidence(esequence,d1endcandt(minid1),d1endcand(minid1),bestFNForwardLengthD1);
        [convectord1,consumd1,flengthd1, solidlengthd1]=computeConfidenceRecursive(esequence,d1endcandt(minid1),d1endcand(minid1),trackingparameters);
        
        v=bbconvector;
        if(size(v,1)>1)
            v=mean(bbconvector);
        end
        
    end
    confidenceData.bestforward1conv=[confidenceData.bestforward1conv;v];
    confidenceData.rforwardd1_conv=[confidenceData.rforwardd1_conv;convectord1];
    confidenceData.rforwardd1_consum=[confidenceData.rforwardd1_consum;consumd1];
    confidenceData.rforwardd1_flength=[confidenceData.rforwardd1_flength;flengthd1];
    confidenceData.rforwardd1_solidlength=[confidenceData.rforwardd1_solidlength;solidlengthd1];
    
end


if(isempty(d2endcand))
    bestFNForwardScoreD2=-1;
    bestFNForwardLengthD2=-1;
    confidenceData.bestforward2conv=[confidenceData.bestforward2conv;zeros(1,6)];
    confidenceData.rforwardd2_conv=[confidenceData.rforwardd2_conv;zeros(1,6)];
    confidenceData.rforwardd2_consum=[confidenceData.rforwardd2_consum;0];
    confidenceData.rforwardd2_flength=[confidenceData.rforwardd2_flength;0];
    confidenceData.rforwardd2_solidlength=[confidenceData.rforwardd2_solidlength;0];
    splitFNMatchScore.forwardd2xy=[splitFNMatchScore.forwardd2xy;-1];
    splitFNMatchScore.forwardd2z=[splitFNMatchScore.forwardd2z;-1];
    splitFNMatchScore.forwardd2gapsize=[splitFNMatchScore.forwardd2gapsize;-1];
    
else
    [scores_gapd2,matchingsd2,playersstartd2,playersendd2]=gapScore(esequence,d2cur,tcur2,d2endcand,d2endcandt,trackingparameters);
    [bestFNForwardScoreD2,minid2]=min(scores_gapd2);
    bestFNForwardLengthD2=d2forwardlengths(minid2);
    
    splitFNMatchScore.forwardd2xy=[splitFNMatchScore.forwardd2xy;distance(esequence{d2endcandt(minid2)}.finalpoints(d2endcand(minid2),1:2)',esequence{tcur2}.finalpoints(d2cur,1:2)')...
        /normv];
    splitFNMatchScore.forwardd2z=[splitFNMatchScore.forwardd2z;abs(esequence{d2endcandt(minid2)}.finalpoints(d2endcand(minid2),3)-esequence{tcur2}.finalpoints(d2cur,3))...
        *trackingparameters.anisotropyvector(3)/normv];
    splitFNMatchScore.forwardd2gapsize=[splitFNMatchScore.forwardd2gapsize;(d2endcandt(minid2)-tcur2)./trackingparameters.interval];
    
    %  'note  forward choice confidences hould be calculated over both branches'
    origparent=esequence{d2endcandt(minid2)}.pred(d2endcand(minid2));
    if(origparent~=-1) %option orig in division 
        parentd1=esequence{d2endcandt(minid2)-1}.suc(origparent,1);
        parentd2=esequence{d2endcandt(minid2)-1}.suc(origparent,2);
        d1lengthforward=traverse_forward(esequence,d2endcandt(minid2),parentd1);
        d2lengthforward=traverse_forward(esequence,d2endcandt(minid2),parentd2);
        targetbiflength=min(d1lengthforward,d2lengthforward); %overwrite length of target with min length of bifurcation
        bestFNForwardLengthD2=targetbiflength;
        [bbconvector,bbconsum]=computeConfidence(esequence,d2endcandt(minid2),parentd1,bestFNForwardLengthD2);
        v=bbconvector;
        [bbconvector,bbconsum]=computeConfidence(esequence,d2endcandt(minid2),parentd2,bestFNForwardLengthD2);
        v=[v;bbconvector];
        if(size(v,1)>1)
            v=mean(v);
        end
        
        [convectord2,consumd2,flengthd2, solidlengthd2]=computeConfidenceRecursive(esequence,d2endcandt(minid2),d2endcand(minid2),trackingparameters);
        
    else %option not orig in division 
        [bbconvector,bbconsum]=computeConfidence(esequence,d2endcandt(minid2),d2endcand(minid2),bestFNForwardLengthD2);
        [convectord2,consumd2,flengthd2, solidlengthd2]=computeConfidenceRecursive(esequence,d2endcandt(minid2),d2endcand(minid2),trackingparameters);   
        v=bbconvector;
        if(size(v,1)>1)
            v=mean(bbconvector); 
        end
    end
    if(size(convectord2,1)>1)
        convectord2=mean(convectord2);
    end
    confidenceData.rforwardd2_conv=[confidenceData.rforwardd2_conv;convectord2];
    confidenceData.rforwardd2_consum=[confidenceData.rforwardd2_consum;consumd2];
    confidenceData.rforwardd2_flength=[confidenceData.rforwardd2_flength;flengthd2];
    confidenceData.rforwardd2_solidlength=[confidenceData.rforwardd2_solidlength;solidlengthd2]; 
    confidenceData.bestforward2conv=[confidenceData.bestforward2conv;v];
    
end
bestEndCandidateInfo(count).t=t;
bestEndCandidateInfo(count).i=i;
bestEndCandidateInfo(count).d1endi=d1cur;
bestEndCandidateInfo(count).d1endt=tcur1;
bestEndCandidateInfo(count).d2endi=d2cur;
bestEndCandidateInfo(count).d2endt=tcur2;
if(isempty(d1endcand))
    bestEndCandidateInfo(count).d1cand=0;
    bestEndCandidateInfo(count).d1candt=0;
else
    bestEndCandidateInfo(count).d1cand=d1endcand(minid1);
    bestEndCandidateInfo(count).d1candt=d1endcandt(minid1);
end
if(isempty(d2endcand))
    bestEndCandidateInfo(count).d2cand=0;
    bestEndCandidateInfo(count).d2candt=0;
else
    bestEndCandidateInfo(count).d2cand=d2endcand(minid2);
    bestEndCandidateInfo(count).d2candt=d2endcandt(minid2);
end


%reset to start for confidence summation etc
d1cur=esequence{t}.suc(i,1);
d2cur=esequence{t}.suc(i,2);
tcur1=esequence{t}.suc_time(i,1);
tcur2=esequence{t}.suc_time(i,2);
d1FP=0;
d2FP=0;


d1conflictend=inf;

for j=1:minlength
    if(isfield(esequence{tcur1},'FP')&isfield(esequence{tcur2},'FP'))
        d1FP=d1FP+esequence{tcur1}.FP(d1cur);
        d2FP=d2FP+esequence{tcur2}.FP(d2cur);
    end
    
    d1curbk=d1cur;
    d2curbk=d2cur;
    d1cur=esequence{tcur1}.suc(d1cur,1);
    d2cur=esequence{tcur2}.suc(d2cur,1);
    tcur1=esequence{tcur1}.suc_time(d1curbk,1);
    tcur2=esequence{tcur2}.suc_time(d2curbk,1);
end

[scores_div,splitscores_div,certainties_div,features]=...
    trackingparameters.DivCostFunction(esequence,i,t,esequence{t}.suc(i,:),esequence{t}.suc_time(i,:), trackingparameters);

BifurcationMeasures=[BifurcationMeasures;features];




    
    % though not needed keep legacy computation of max
    % length available option
    if(~isempty(FNbackcand1lengths))
        FNbackcand1lengths=max(FNbackcand1lengths);
    else
        FNbackcand1lengths=-1;
    end
    if(~isempty(FNbackcand2lengths))
        FNbackcand2lengths=max(FNbackcand2lengths);
    else
        FNbackcand2lengths=-1;
    end

if trackingparameters.recordanswers
    
       %obsolete variables that are part of answer vector
    FNbackmutual=-1;
    d2conflictend=inf;d1conflictend=inf;
    d1FPcon=0;d2FPcon=0;d1consumconfidence=0;d2consumconfidence=0;d3FPcon=0;d3consumconfidence=0;
    conflictendlength=-1;
    d1sumconfidence=0;
    d2sumconfidence=0;
    
    
    %correct
    canswer_suc=esequence{t}.correct_suc(i,:);
    canswer_suc_time=esequence{t}.correct_suc_time(i,:);
    divanswerright=min((canswer_suc==esequence{t}.suc(i,:)))&min(canswer_suc_time==esequence{t}.suc_time(i));
    flip=esequence{t}.suc(i,:);
    flip(1)=esequence{t}.suc(i,2); flip(2)=esequence{t}.suc(i,1);
    flipt(1)=esequence{t}.suc_time(i,2); flipt(2)=esequence{t}.suc_time(i,1);
    divanswerright=divanswerright|(min((canswer_suc==flip))&min(canswer_suc_time==flipt));
    
    d1length=d1length./trackingparameters.interval;
    d2length=d2length./trackingparameters.interval;
    
    if(bestdaughter~=-1)
        if (bestdaughter==1)
            [match,matchi,matchj ] = checkIfRealSuccessors(d1candt(bestIndex),d1cand(bestIndex),esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),esequence );
        else
            [match,matchi,matchj ] = checkIfRealSuccessors(d2candt(bestIndex),d2cand(bestIndex),esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),esequence );
            
        end
        simpleFNcorrect=[simpleFNcorrect;match];
    else
        simpleFNcorrect=[simpleFNcorrect;-1];
    end
    
 

    removed=[removed;t,i,d1length,d2length,splitscores_div(1:3),divanswerright,...
        d1sumconfidence,d2sumconfidence,d1FP,d2FP,d1conflictend,d2conflictend,d1consumconfidence,...
        d2consumconfidence,d1FPcon,d2FPcon,d3consumconfidence,d3FPcon, conflictendlength,dFNf,dFNb,...
        FNbackcand1,FNbackcand2,branchendFNoption1,branchendFNoption2,FNbackcand1lengths,FNbackcand2lengths,...
        bestFNBackScore,bestFNBackLength,bestFNForwardScoreD1,bestFNForwardScoreD2,...
        bestFNForwardLengthD1,bestFNForwardLengthD2, bestFNBackCorrect,FNbackmutual];
end

count=count+1;
%coallate data for indivdual classifications
minsize=min(d1length,d2length);
%collect data for whichever daughter is smaller
if(minsize==d1length)
    minsizedaughterbestscorelength=bestFNForwardLengthD1;
    minsizedaughterbestscore=bestFNForwardScoreD1;
    minbranchforwardconv=confidenceData.bestforward1conv(count-1,:);
    minbranchforwardconv_recursive=confidenceData.rforwardd1_conv(count-1,:);
    minbranchforwardhops_recursive=confidenceData.rforwardd1_consum(count-1,:);
    minbranchforward_solidlength=confidenceData.rforwardd1_solidlength(count-1);
    minbranchforward_timelength=confidenceData.rforwardd1_flength(count-1);
    minbranchforwardxy=splitFNMatchScore.forwardd1xy(count-1);
    minbranchforwardz=splitFNMatchScore.forwardd1z(count-1);
    minbranchforwardgapsize=splitFNMatchScore.forwardd1gapsize(count-1);
else
    minsizedaughterbestscorelength=bestFNForwardLengthD2;
    minsizedaughterbestscore=bestFNForwardScoreD2;
    minbranchforwardconv=confidenceData.bestforward2conv(count-1,:);
    minbranchforwardconv_recursive=confidenceData.rforwardd2_conv(count-1,:);
    minbranchforwardhops_recursive=confidenceData.rforwardd2_consum(count-1,:);
    minbranchforward_solidlength=confidenceData.rforwardd2_solidlength(count-1);
    minbranchforward_timelength=confidenceData.rforwardd2_flength(count-1);
    minbranchforwardxy=splitFNMatchScore.forwardd2xy(count-1);
    minbranchforwardz=splitFNMatchScore.forwardd2z(count-1);
    minbranchforwardgapsize=splitFNMatchScore.forwardd2gapsize(count-1);
end

ncells=mean(esequence{t}.selfdistance)./mean(esequence{t}.finaldiams);

alldaughterdata=[BifurcationMeasures(count-1,:),confidenceData.bifconvector(count-1,:),minsize,ncells];
allforwarddata=[minbranchforwardgapsize,minsizedaughterbestscore',minsizedaughterbestscorelength',minbranchforwardconv,minbranchforwardxy,minbranchforwardz,minbranchforward_solidlength,minbranchforward_timelength];
allbackdata=[splitFNMatchScore.backgapsize(count-1,1),bestFNBackScore,bestFNBackLength,confidenceData.bestbackconv(count-1,:),splitFNMatchScore.backxy(count-1,:),splitFNMatchScore.backz(count-1)];
end

