function  [ data]= computeNonDivisionConfidenceVector(esequence, t, i,trackingparameters)
%computes feature vector for the confidence of a bifurcation
%assumes the passed in nucleus is about to bifurcate
data=calculateCellPairVectorNondivision_wdiam(esequence{t},i,esequence{esequence{t}.suc_time(i,1)},esequence{t}.suc(i,1),trackingparameters.anisotropyvector);

%wide window measures
data=[data,calculateWideWindowCellPairVectorNondivision(esequence,t,i,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),trackingparameters.anisotropyvector,trackingparameters.wideWindow)];

%alternatives check
suc=esequence{t}.suc(i,1);
suc_t=esequence{t}.suc_time(i,1);


%I took filter out because it made results worse, though should probably
%recheck
%data=[data,filterCandidateLength(esequence{t}.sucessor_suitors{i},esequence{t+1}.delete)];
%data=[data,filterCandidateLength(esequence{suc_t}.predecessor_suitors{suc},esequence{suc_t-1}.delete)];

%4 counts of nn claimants at 2 endpoints
data=[data,length(esequence{t}.sucessor_suitors{i})];
data=[data,length(esequence{suc_t}.predecessor_suitors{suc})];
data=[data,length(esequence{t}.predecessor_suitors{i})];
data=[data,length(esequence{suc_t}.sucessor_suitors{suc})];

% 2 interior comparisons btw nn and link distance
% data=[data,nnInterfeerenceConfidenceFeature( esequence,t,i,0,trackingparameters )];
%data=[data,nnInterfeerenceConfidenceFeature( esequence,t,i,1,trackingparameters )];
 
%4 distance to endpoint tests at endpoint times 
data=[data,endInterfeerenceConfidenceFeature(esequence, t,t,i,0,trackingparameters )];
data=[data,endInterfeerenceConfidenceFeature(esequence, suc_t,suc_t,suc,1,trackingparameters )];

data=[data,endInterfeerenceConfidenceFeature(esequence, t,t,i,1,trackingparameters )];
data=[data,endInterfeerenceConfidenceFeature(esequence, suc_t,suc_t,suc,0,trackingparameters )];
%{

%forward
data=[data,endInterfeerenceConfidenceFeature(esequence, t+1,t,i,0,trackingparameters )];
data=[data,endInterfeerenceConfidenceFeature(esequence, suc_t+1,suc_t,suc,1,trackingparameters )];

data=[data,endInterfeerenceConfidenceFeature(esequence, t+1,t,i,1,trackingparameters )];
data=[data,endInterfeerenceConfidenceFeature(esequence, suc_t+1,suc_t,suc,0,trackingparameters )];
%back
data=[data,endInterfeerenceConfidenceFeature(esequence, t-1,t,i,0,trackingparameters )];
data=[data,endInterfeerenceConfidenceFeature(esequence, suc_t-1,suc_t,suc,1,trackingparameters )];

data=[data,endInterfeerenceConfidenceFeature(esequence, t-1,t,i,1,trackingparameters )];
data=[data,endInterfeerenceConfidenceFeature(esequence, suc_t-1,suc_t,suc,0,trackingparameters )];

%}
end

