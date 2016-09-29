function [scores,scores2,scores3,scores4]=distanceCostFunction(esequence,i,t,candidates,candidate_times,trackingparameters)
%Score potential matches based on product of gfp and distance
%scores2 are split scores 3 are dummy to match model based function that
%returns confidence
scores=zeros(size(candidates));
c1=zeros(size(candidates));
c2=zeros(size(candidates));
for j=1:length(candidates)
    scores(j)=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidate_times(j)}.finalpoints(candidates(j),:)',trackingparameters.anisotropyvector);
    if(isfield(trackingparameters,'abscutoff')&&trackingparameters.abscutoff)
    else
    scores(j)=scores(j)./mean(esequence{t}.selfdistance);
    end
    
    c1(j)=esequence{t}.confidences(i);
     c2(j)=esequence{candidate_times(j)}.confidences(candidates(j));


end

scores2=[scores,c1,c2];
scores3=scores;
scores4=scores;
end

