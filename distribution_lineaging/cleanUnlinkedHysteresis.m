function esequence=cleanUnlinkedHysteresis(esequence,trackingparameters);
%walks through cells and marks as dead all those unlinked and below high
%threshold
'cleaning unlinked dim'
for t=trackingparameters.starttime:trackingparameters.endtime
    for i=1:size(esequence{t}.finalpoints,1)
        
        hysteresisbad=...
            esequence{t}.finalmaximas(i)<trackingparameters.hysteresis_intensityhigh&&...
            esequence{t}.pred(i)==-1&&esequence{t}.suc(i,1)==-1;
        
        
        if(hysteresisbad)
            esequence{t}.delete(i)=1;
        end
    end
end


end

