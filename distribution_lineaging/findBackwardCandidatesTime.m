function [candidates,candidates_t]=findBackwardCandidatesTime(esequence,i,t,trackingparameters);
      
% finds possible backward gap matches of start i at time t at
%previous ends
candidates=[];
candidates_t=[];
temporalcutoff=trackingparameters.temporalcutoff;
%start of 2 implies no zero length FN  1 implies allowed
%'allowing zero length fn'

for offset=trackingparameters.temporalcutoffstart:temporalcutoff
    if(t-offset>=1&&~isempty(esequence{t-offset}.finalpoints))
       
        spatialcutoff=trackingparameters.forwardcutoff(t-offset);
        distances=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{t-offset}.finalpoints',trackingparameters.anisotropyvector);
        %within cutoff and a midpoint
        ccandidates=find(~esequence{t-offset}.delete&...
            distances'<spatialcutoff&esequence{t-offset}.suc(:,1)==-1);
        candidates=[candidates;ccandidates];
        candidates_t=[candidates_t;ones(size(ccandidates))*t-offset];
    end
end

end

