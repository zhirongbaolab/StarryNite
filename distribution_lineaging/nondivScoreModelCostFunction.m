function [scores,splitscores,certainties_nondiv,nondivscorecomponents]=nondivScoreModelCostFunction(esequence,i,t,candidates,candidate_times,trackingparameters)
%Score potential division matches based on product of scores for triple and pairs

scores=zeros(size(candidates,1),1);
for j=1:size(candidates,1)
    
      pairdata=calculateCellPairVectorNondivision(esequence{t},i,esequence{candidate_times(j)},candidates(j),trackingparameters.anisotropyvector);
    score=(1./mvnpdf(pairdata,trackingparameters.model.nodiv_mean,trackingparameters.model.nodiv_std));
    scores(j)=log(score);
     
end
splitscores=scores;
certainties_nondiv=scores;
nondivscorecomponents=scores;
end

