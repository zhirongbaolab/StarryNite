function [scores,splitscores,certainties,features]=divScoreModelCostFunction(esequence,i,t,candidates,candidate_times,trackingparameters)
%Score potential division matches based on product of scores for triple and pairs
features=[];
scores=zeros(size(candidates,1),1);
certainties=zeros(size(candidates,1),1);
splitscores=zeros(size(candidates,1),6);
for j=1:size(candidates,1)
    tripledata=calculateCellTripleVector(esequence{t},i,esequence{candidate_times(j,1)},candidates(j,1),esequence{candidate_times(j,2)},candidates(j,2),trackingparameters.anisotropyvector,trackingparameters.interval);
    agree=1./mvnpdf(tripledata,trackingparameters.model.div_triple_mean,trackingparameters.model.div_triple_std);
    
      
    pairdata1=calculateCellPairVector(esequence{t},i,esequence{candidate_times(j,1)},candidates(j,1),trackingparameters.anisotropyvector);
    score1=(1./mvnpdf(pairdata1,trackingparameters.model.div_mean,trackingparameters.model.div_std));
    
    pairdata2=calculateCellPairVector(esequence{t},i,esequence{candidate_times(j,2)},candidates(j,2),trackingparameters.anisotropyvector);
    score2=(1./mvnpdf(pairdata2,trackingparameters.model.div_mean,trackingparameters.model.div_std));
    splitscores(j,1)=agree;
    splitscores(j,2)=score1;
    splitscores(j,3)=score2;
    
    scores(j)=log(agree*score1*score2);
      features=[features;pairdata1,pairdata2,tripledata];

end
end

