function [candidates,candidates_t]=findBackwardCandidatesTime_NN(esequence,i,t,trackingparameters);
      
%findBackwardCandidates finds possible backward matches of start i at time t at
%previous t-1 via x nn backward
candidates=[];
candidates_t=[];
temporalcutoff=trackingparameters.temporalcutoff;
for offset=1:temporalcutoff
    if(t-offset>=1&&~isempty(esequence{t-offset}.finalpoints))
        spatialcutoff=trackingparameters.forwardcutoff(t-offset);
        %spatialcutoff=inf
        distances=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{t-offset}.finalpoints',trackingparameters.anisotropyvector);
        ccandidates=[];
        number=trackingparameters.nnnumber;
        if(offset>1)
            number=trackingparameters.nnnumber_gap;
        end
        for j=1:min(number,length(distances))
            [v,im]=min(distances);
            if(v<spatialcutoff&&esequence{t-offset}.suc(im,2)==-1) %if there is someplace to put it and  is within thresh      
                % if((j<2||v<spatialcutoff)&esequence{t-offset}.suc(im,2)==-1) %if there is someplace to put it and  is within thresh
                ccandidates=[ccandidates;im];
            end
            distances(im)=Inf;
        end
        candidates=[candidates;ccandidates];
        candidates_t=[candidates_t;ones(size(ccandidates))*t-offset];
    end
end

end

