function [scores,splitscores,confidences]=divDistanceCostFunction(esequence,i,t,candidates,candidate_times,trackingparameters)
%Score potential matches based on product of gfp and distance

scores=zeros(size(candidates,1),1);
for j=1:size(candidates,1)
    midpoint=esequence{candidate_times(j,1)}.finalpoints(candidates(j,1),:)+esequence{candidate_times(j,2)}.finalpoints(candidates(j,2),:);
    midpoint=midpoint./2;
    %dis1=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidate_times(j,1)}.finalpoints(candidates(j,1),:)',trackingparameters.anisotropy);
    %dis2=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidate_times(j,2)}.finalpoints(candidates(j,2),:)',trackingparameters.anisotropy);
    
   % gfpratio=esequence{t}.totalGFP(i)/(esequence{candidate_times(j,1)}.totalGFP(candidates(j,1))+esequence{candidate_times(j,2)}.totalGFP(candidates(j,2)));
   % gfpratio=max(gfpratio,1/gfpratio);
   %scores(j)=(dis1+dis2);%*gfpratio;
   scores(j)=distance_anisotropic(esequence{t}.finalpoints(i,:)',midpoint',trackingparameters.anisotropyvector);
 % scores(j)=scores(j)./mean(esequence{t}.finaldiams);
 %mode in which distance normalization is ignored   
 if(isfield(trackingparameters,'abscutoff')&&trackingparameters.abscutoff)
     else
         scores(j)=scores(j)./mean(esequence{t}.selfdistance);
     end
end
confidences=scores;
splitscores=[scores,scores,scores,scores,scores,scores];%dummy
features=scores;
end

