function [esequence,numplayers,endplayerssameasdiv ]= processFNBifurcation...
    ( esequence,t,i,bestmatchings, bestmatchingsplayerstart,...
    bestmatchingplayerend, bestdaughter,bestIndex,d1cand,d2cand,...
    d1candt,d2candt,trackingparameters)
%process a FN bifurcation by linknig if 1 end to 1 end dirty case or
%relinknig all conflicting players in best configuration given by
%bestmatchings if clean conflict case (2 or 3 players)

matchinfo=bestmatchings;
playersstartinfo=bestmatchingsplayerstart;
playersendinfo=bestmatchingplayerend;
finaltime=t+1;
if(bestdaughter==1)
    starttime=d1candt(bestIndex);
else
    starttime=d2candt(bestIndex);
end

% ['start of FN link ',num2str(matchinfo)]
numplayers=length(find(playersendinfo));

%check whether the two conflict candidates are identical to the division
%if not this is not a clean case strictly speaking 
d1ispresent=max(esequence{t}.suc(i,1)==playersendinfo(1:2));
d2ispresent=max(esequence{t}.suc(i,2)==playersendinfo(1:2));
endplayerssameasdiv=d1ispresent&d2ispresent;

%similarly check that one of the conflic start players is actually
%connected to something otherwise it is also not a strict clean case
if (numplayers>1)
    startlinkpresent=(max(esequence{starttime}.suc(playersstartinfo(1:2),1))>-1);
end

if(bestdaughter==1)
    loosestart=d1cand(bestIndex);
else
    loosestart=d2cand(bestIndex);
end

%if not clean in matching # at conflict sense or in sense of the conflict
%not overlapping with the existing links, treat it as a simple 1-1 FN link
if (numplayers==1||~endplayerssameasdiv||~startlinkpresent)
    if(bestdaughter==1)
        %{
        FNmade=[FNmade; d1candt(bestIndex),esequence{t}.suc_time(i,1)...
            esequence{d1candt(bestIndex)}.finalpoints(d1cand(bestIndex),:),...
            esequence{esequence{t}.suc_time(i,1)}.finalpoints(esequence{t}.suc(i,1),:)];
        %}
        esequence{d1candt(bestIndex)}.suc(d1cand(bestIndex),1)=esequence{t}.suc(i,1);
        esequence{d1candt(bestIndex)}.suc_time(d1cand(bestIndex),1)=esequence{t}.suc_time(i,1);
        esequence{esequence{t}.suc_time(i,1)}.pred(esequence{t}.suc(i,1))=d1cand(bestIndex);
        esequence{esequence{t}.suc_time(i,1)}.pred_time(esequence{t}.suc(i,1))=d1candt(bestIndex);
        esequence{t}.suc_time(i,1)=esequence{t}.suc_time(i,2);
        esequence{t}.suc(i,1)=esequence{t}.suc(i,2);
        
        esequence{t}.suc_time(i,2)=-1;
        esequence{t}.suc(i,2)=-1;
    else
        %{
        FNmade=[FNmade; d2candt(bestIndex),esequence{t}.suc_time(i,2)...
            esequence{d2candt(bestIndex)}.finalpoints(d2cand(bestIndex),:),...
            esequence{esequence{t}.suc_time(i,2)}.finalpoints(esequence{t}.suc(i,2),:)];
        %}
        esequence{d2candt(bestIndex)}.suc(d2cand(bestIndex),1)=esequence{t}.suc(i,2);
        esequence{d2candt(bestIndex)}.suc_time(d2cand(bestIndex),1)=esequence{t}.suc_time(i,2);
        esequence{esequence{t}.suc_time(i,2)}.pred(esequence{t}.suc(i,2))=d2cand(bestIndex);
        esequence{esequence{t}.suc_time(i,2)}.pred_time(esequence{t}.suc(i,2))=d2candt(bestIndex);
        
        esequence{t}.suc_time(i,2)=-1;
        esequence{t}.suc(i,2)=-1;
    end
else
    if( numplayers==2)
        'clean 2 case';
        oneback=esequence{finaltime}.pred(playersendinfo(1));
        % non -1 start player suc
        oneforward=max(esequence{starttime}.suc(playersstartinfo(1:2)));
        distancepair1=distance_anisotropic(esequence{starttime}.finalpoints(playersstartinfo(1),:)',...
            esequence{starttime+1}.finalpoints(oneforward,:)',trackingparameters.anisotropyvector)+...
            distance_anisotropic(esequence{finaltime}.finalpoints(playersendinfo(matchinfo(1)),:)',...
            esequence{finaltime-1}.finalpoints(oneback,:)',trackingparameters.anisotropyvector);                                        ;
        distancepair2=distance_anisotropic(esequence{starttime}.finalpoints(playersstartinfo(2),:)',...
            esequence{starttime+1}.finalpoints(oneforward,:)',trackingparameters.anisotropyvector)+...
            distance_anisotropic(esequence{finaltime}.finalpoints(playersendinfo(matchinfo(2)),:)',...
            esequence{finaltime-1}.finalpoints(oneback,:)',trackingparameters.anisotropyvector);
        if(distancepair1<distancepair2)
            gapstart=playersstartinfo(2);
            gapend=playersendinfo(matchinfo(2));
            %link pair 1 to middle
            startmiddleclaimant=playersstartinfo(1);
            endmiddleclaimant=playersendinfo(matchinfo(1));
        else
            gapstart=playersstartinfo(1);
            gapend=playersendinfo(matchinfo(1));
            %link pair 2 to middle
            startmiddleclaimant=playersstartinfo(2);
            endmiddleclaimant=playersendinfo(matchinfo(2));
        end
        %link better claimants to middle
        esequence{starttime}.suc(startmiddleclaimant,1)=oneforward;
        esequence{starttime}.suc_time(startmiddleclaimant,1)=starttime+1;
        esequence{starttime+1}.pred(oneforward)=startmiddleclaimant;
        esequence{starttime+1}.pred_time(oneforward)=starttime;
        esequence{finaltime-1}.suc(oneback,1)= endmiddleclaimant;
        esequence{finaltime-1}.suc_time(oneback,1)=finaltime;
        esequence{finaltime}.pred(endmiddleclaimant)=oneback;
        esequence{finaltime}.pred_time(endmiddleclaimant)=finaltime-1;
        %link gap pair to eachother
        esequence{starttime}.suc(gapstart,1)=gapend;
        esequence{starttime}.suc_time(gapstart,1)=finaltime;
        esequence{finaltime}.pred(gapend)=gapstart;
        esequence{finaltime}.pred_time(gapend)=starttime;
        esequence{finaltime-1}.suc(oneback,2)= -1;
        esequence{finaltime-1}.suc_time(oneback,2)=-1;
        
    else
        if(numplayers==3)
            
            'clean 3 case';
            matchofstart3=playersendinfo(matchinfo(3));
            if(matchofstart3==esequence{t}.suc(i,1))
                %d1 is match
                %d2 is gap
                gapend=esequence{t}.suc(i,2);
                %cut off d2
                esequence{t}.suc_time(i,2)=-1;
                esequence{t}.suc(i,2)=-1;
                
            else
                %d2 is match
                %d1 is gap
                gapend=esequence{t}.suc(i,1);
                esequence{t}.suc_time(i,1)=esequence{t}.suc_time(i,2);
                esequence{t}.suc(i,1)=esequence{t}.suc(i,2);
                esequence{t}.suc_time(i,2)=-1;
                esequence{t}.suc(i,2)=-1;
            end
            
            
            matchofend3i=find(matchinfo==3);
            if (matchofend3i==3)
                'warning: clean 3 FN conflict was matched as breaking nonconflicting link arbitrary solution being executed'
            end
            if(matchofend3i==1)
                othermatch=playersstartinfo(2);
            else
                othermatch=playersstartinfo(1);
            end
            matchofend3=playersstartinfo(find(matchinfo==3));
            if(matchofend3==loosestart)
              
                
                %loose is match swap loose and
                %other so that other is loose
                gapstart=othermatch;
                %real match as sucessor in chain
                esequence{starttime}.suc(loosestart,1)=esequence{starttime}.suc(othermatch,1);
                esequence{starttime}.suc_time(loosestart,1)=esequence{starttime}.suc_time(othermatch,1);
                %point chain at real successor
                esequence{esequence{starttime}.suc_time(loosestart,1)}.pred(esequence{starttime}.suc(loosestart,1))=loosestart;
                esequence{esequence{starttime}.suc_time(loosestart,1)}.pred_time(esequence{starttime}.suc(loosestart,1))=starttime;
                %&cut other end loose
                esequence{starttime}.suc(othermatch,1)=-1;
                esequence{starttime}.suc_time(othermatch,1)=-1;
                
            else
                gapstart=loosestart;
            end
            %link start and end which are now
            %both guaranteed to be loose
            esequence{starttime}.suc(gapstart,1)=gapend;
            esequence{starttime}.suc_time(gapstart,1)=finaltime;
            esequence{finaltime}.pred(gapend)=gapstart;
            esequence{finaltime}.pred_time(gapend)=starttime;
            
        else
            'warning: fn class is not any of known topology types, shouldnt happen' 
        end
    end
end
end

