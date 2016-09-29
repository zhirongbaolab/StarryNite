function  [ datavector]= computeDeathConfidenceVector(esequence, t, i,trackingparameters)
%computes feature vector for the confidence of a bifurcation
%assumes the passed in nucleus is about to bifurcate
datavector=[esequence{t}.totalGFP(i)/mean(esequence{t}.totalGFP),...
    esequence{t}.avgGFP(i)/mean(esequence{t}.avgGFP),...
    esequence{t}.finaldiams(i)/mean(esequence{t}.finaldiams)];

[depth ] = traverse_backdivstop( esequence,t,i );
tp=depth;
if(depth>10)
    tp=tp-9;
end
ic=i;
tc=t;
for j=0:tp-2
    ict=ic;
    ic=esequence{tc}.pred(ic);
    tc=esequence{tc}.pred_time(ict);
end

datavectorearly=[esequence{tc}.totalGFP(ic)/mean(esequence{tc}.totalGFP),...
    esequence{tc}.avgGFP(ic)/mean(esequence{tc}.avgGFP),...
    esequence{tc}.finaldiams(ic)/mean(esequence{tc}.finaldiams)];

datavector=[datavector,datavectorearly];

datavector=[datavector,endInterfeerenceConfidenceFeature(esequence, t,t,i,0,trackingparameters )];
datavector=[datavector,endInterfeerenceConfidenceFeature(esequence, t,t,i,1,trackingparameters )];

%{
datavector=[datavector,endInterfeerenceConfidenceFeature(esequence, t+1,t,i,0,trackingparameters )];
datavector=[datavector,endInterfeerenceConfidenceFeature(esequence, t+1,t,i,1,trackingparameters )];
datavector=[datavector,endInterfeerenceConfidenceFeature(esequence, t-1,t,i,0,trackingparameters )];
datavector=[datavector,endInterfeerenceConfidenceFeature(esequence, t-1,t,i,1,trackingparameters )];

%}
datavector(isinf(datavector))=0;
datavector(isnan(datavector))=0;
end

