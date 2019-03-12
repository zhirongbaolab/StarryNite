function [candidates,candidates_t]=findForwardCandidatesTime(esequence,i,t,trackingparameters)

%findForwardCandidates finds possible forward matches of i at time t at
%future ts
candidates=[];
candidates_t=[];
temporalcutoff=trackingparameters.temporalcutoff;
for offset=2:min(temporalcutoff,trackingparameters.endtime-t)
    if(~isempty(esequence{t+offset}.finalpoints))
        
        spatialcutoff=trackingparameters.forwardcutoff(t+offset);
  
            distances=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{t+offset}.finalpoints',trackingparameters.anisotropyvector);
  
        %within cutoff and a potential start in that its predecessor has 2
        %successors
        preddiv=zeros(size(esequence{t+offset}.pred));
        
        %has a pred and it is at previous timepoint
        haspred=esequence{t+offset}.pred~=-1&esequence{t+offset}.pred_time==t+offset-1;
        %set those that havepred to true if their predecessor has a second
        %sucessor
        preddiv(haspred)=esequence{t+offset-1}.suc(esequence{t+offset}.pred(haspred),2)~=-1;
        
        %previous inefficient implementation should be identical given fact
        %that divisions never cover gaps in current implementation
        %for j=1:length(preddiv)
        %    if(esequence{t+offset}.pred(j)~=-1)
        %   preddiv(j)=esequence{esequence{t+offset}.pred_time(j)}.suc(esequence{t+offset}.pred(j),2)~=-1;
        %    end
        %end
        ccandidates=find((~esequence{t+offset}.delete')&(distances<spatialcutoff)&(esequence{t+offset}.pred'==-1|preddiv'));
        candidates=[candidates;ccandidates'];
        candidates_t=[candidates_t;ones(size(ccandidates'))*t+offset];
    end
end

end

