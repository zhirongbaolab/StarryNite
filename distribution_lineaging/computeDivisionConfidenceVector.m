function  [ divdata]= computeDivisionConfidenceVector(esequence, t, i,trackingparameters)
%computes feature vector for the confidence of a bifurcation
%assumes the passed in nucleus is about to bifurcate
pairdata=calculateCellPairVector(esequence{t},i,esequence{esequence{t}.suc_time(i,1)},esequence{t}.suc(i,1),trackingparameters.anisotropyvector);
pairdata2=calculateCellPairVector(esequence{t},i,esequence{esequence{t}.suc_time(i,2)},esequence{t}.suc(i,2),trackingparameters.anisotropyvector);
tripledata=calculateCellTripleVector(esequence{t},i...
    ,esequence{esequence{t}.suc_time(i,1)},esequence{t}.suc(i,1),...
    esequence{esequence{t}.suc_time(i,2)},esequence{t}.suc(i,2),...
    trackingparameters.anisotropyvector,trackingparameters.interval);

divdata=[pairdata,pairdata2,tripledata];

divdata=[divdata,calculateWideWindowCellTripleVector(...
    esequence,t,i,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),...
    esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),...
    trackingparameters.anisotropyvector,trackingparameters.wideWindow)];


%alternatives check
suc=esequence{t}.suc(i,1);
suc2=esequence{t}.suc(i,2);
suc_t=esequence{t}.suc_time(i,1);
suc2_t=esequence{t}.suc_time(i,2);
%divdata=[divdata,filterCandidateLength(esequence{t}.sucessor_suitors{i},esequence{t+1}.delete)];
%divdata=[divdata,filterCandidateLength(esequence{suc_t}.predecessor_suitors{suc},esequence{suc_t-1}.delete)];
%divdata=[divdata,filterCandidateLength(esequence{suc2_t}.predecessor_suitors{suc2},esequence{suc2_t}.delete)];


divdata=[divdata,length(esequence{t}.sucessor_suitors{i})];
divdata=[divdata,length(esequence{suc_t}.predecessor_suitors{suc})];
divdata=[divdata,length(esequence{suc2_t}.predecessor_suitors{suc2})];

divdata=[divdata,length(esequence{t}.predecessor_suitors{i})];
divdata=[divdata,length(esequence{suc_t}.sucessor_suitors{suc})];
divdata=[divdata,length(esequence{suc2_t}.sucessor_suitors{suc2})];

%new calc apandoned
%{
% additional nn claimant vs link computes min for both daughters in case of
% division
 divdata=[divdata,nnInterfeerenceConfidenceFeature( esequence,t,i,0,trackingparameters )];
 divdata=[divdata,nnInterfeerenceConfidenceFeature( esequence,t,i,1,trackingparameters )];
 %}

%initial
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, t,t,i,0,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc_t,suc_t,suc,1 ,trackingparameters)];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc2_t,suc2_t,suc2,1 ,trackingparameters)];
%self
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, t,t,i,1,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc_t,suc_t,suc,0,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc2_t,suc2_t,suc2,0,trackingparameters )];

%{
%forward in time
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, t+1,t,i,0,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc_t+1,suc_t,suc,1 ,trackingparameters)];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc2_t+1,suc2_t,suc2,1 ,trackingparameters)];
%self
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, t+1,t,i,1,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc_t+1,suc_t,suc,0,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc2_t+1,suc2_t,suc2,0,trackingparameters )];

%back in time
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, t-1,t,i,0,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc_t-1,suc_t,suc,1 ,trackingparameters)];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc2_t-1,suc2_t,suc2,1 ,trackingparameters)];
%self
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, t-1,t,i,1,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc_t-1,suc_t,suc,0,trackingparameters )];
divdata=[divdata,endInterfeerenceConfidenceFeature(esequence, suc2_t-1,suc2_t,suc2,0,trackingparameters )];

%}
end

