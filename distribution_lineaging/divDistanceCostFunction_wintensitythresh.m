function [scores,splitscores,confidences]=divDistanceCostFunction_wintensitythresh(esequence,i,t,candidates,candidate_times,trackingparameters)
%Score potential matches based on product of gfp and distance

scores=zeros(size(candidates,1),1);
for j=1:size(candidates,1)
    midpoint=esequence{candidate_times(j,1)}.finalpoints(candidates(j,1),:)+esequence{candidate_times(j,2)}.finalpoints(candidates(j,2),:);
    midpoint=midpoint./2;
    %dis1=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidate_times(j,1)}.finalpoints(candidates(j,1),:)',trackingparameters.anisotropy);
    %dis2=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidate_times(j,2)}.finalpoints(candidates(j,2),:)',trackingparameters.anisotropy);
    
   gfpratio1=esequence{candidate_times(j,1)}.totalGFP(candidates(j,1))/esequence{t}.totalGFP(i);
   gfpratio2=esequence{candidate_times(j,2)}.totalGFP(candidates(j,1))/esequence{t}.totalGFP(i);

        
   % gfpratio=max(gfpratio,1/gfpratio);
   %scores(j)=(dis1+dis2);%*gfpratio;
   scores(j)=distance_anisotropic(esequence{t}.finalpoints(i,:)',midpoint',trackingparameters.anisotropyvector);
if(gfpratio1>.775||gfpratio1<.35||gfpratio2>.775||gfpratio2<.35)
  scores(j)=inf;
end
   % scores(j)=scores(j)./mean(esequence{t}.finaldiams);
 %mode in which distance normalization is ignored   
 if(isfield(trackingparameters,'abscutoff')&&trackingparameters.abscutoff)
     else
         scores(j)=scores(j)./mean(esequence{t}.selfdistance);
 end
     

end
confidences=scores;
splitscores=[scores,scores,scores,scores,scores,scores];%dummy
end

