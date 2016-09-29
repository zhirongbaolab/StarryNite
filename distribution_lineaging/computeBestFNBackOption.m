function [bestdaughter,bestIndex,bestFNBackScore,iscorrect,bestmatchings,bestmatchingsplayerstart,bestmatchingplayerend]=...
    computeBestFNBackOption(d1cand,d1candt,d2cand,d2candt,tcur1,d1cur,d2cur,esequence,trackingparameters);
%return the best FN backward option for the division with d1 and d2, the candidates have already been extracted
%each needs to be pulled out and scored with the min score over both
%daughters returned along with the ID of the daughter and the index into
%the options returns -1's if the candidate options are empty
allscores=[];
allindicies=[];
allt=[];
allplayersstart=[];

allplayersend=[];
allmatchings=[];
%given that scoring is one to many call once for each backward candidate
%with origin
for candidate=1:length(d1cand)
    [scores_gap,matchings,playersstart,playersend]=gapScore(esequence,d1cand(candidate),d1candt(candidate),d1cur,tcur1,trackingparameters);
    allscores=[allscores;scores_gap];
    allindicies=[allindicies;1,candidate];
    allt=[allt;d1candt(candidate)];
    allplayersstart=[allplayersstart;playersstart];
    
    allplayersend=[allplayersend;playersend];
    
    allmatchings=[allmatchings;matchings];
end

for candidate=1:length(d2cand)
    [scores_gap,matchings,playersstart,playersend]=gapScore(esequence,d2cand(candidate),d2candt(candidate),d2cur,tcur1,trackingparameters);
    allscores=[allscores;scores_gap];
    allindicies=[allindicies;2,candidate];
    allt=[allt;d2candt(candidate)];
    allplayersstart=[allplayersstart;playersstart];
    
    allplayersend=[allplayersend;playersend];
    
    allmatchings=[allmatchings;matchings];
end
if(isempty(allscores))
    bestdaughter=-1;
    bestIndex=-1;
    bestFNBackScore=-1;
    iscorrect=0;
    bestmatchings=-1;
    bestmatchingsplayerstart=-1;
    bestmatchingplayerend=-1;
else
    [bestFNBackScore,i]=min(allscores);
    bestdaughter=allindicies(i,1);
    bestIndex=allindicies(i,2);
    bestmatchings=allmatchings(i,:);
    bestmatchingsplayerstart=allplayersstart(i,:);
    bestmatchingplayerend=allplayersend(i,:);
    
    gapchoicescorrect=zeros(1,3);
    if trackingparameters.recordanswers
    for k=1:3
        if(allplayersstart(i,k)~=0)
            [match,matchi,matchj ] = checkIfRealSuccessors(allt(i),allplayersstart(i,k),tcur1,allplayersend(i,allmatchings(i,k)),esequence );
        else
            match=1;
        end
        gapchoicescorrect(1,k)=match;
    end
    end
    iscorrect=sum(gapchoicescorrect)==3;
    
end



end

