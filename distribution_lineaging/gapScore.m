function [scores_gap,matchings,playersstart,playersend]=gapScore(esequence,i,t,candidates,candidates_t,trackingparameters);
% given canididate free ends finds conflicting other players
%scores based no considering 1-3 permutations of FN/normal linking among
%players
averagediam=mean(esequence{t}.selfdistance);
scores_gap=zeros(size(candidates,1),1);
matchings=zeros(size(candidates,1),3);
playersstart=zeros(size(candidates,1),3);
playersend=zeros(size(candidates,1),3);
for j=1:size(candidates_t,1)
    [ startconflictplayers,startbacktrace,endconflictplayers,endforwardtrace ] = ...
        FindFNplayers( t,i,candidates_t(j),candidates(j),esequence );
    clean=(length(startconflictplayers)==2&length(startconflictplayers)==length(endconflictplayers))&...
        length(endforwardtrace)==length(startbacktrace)&...
        (isempty(startbacktrace)||startbacktrace~=-1)&(isempty(endforwardtrace)||endforwardtrace~=-1);

    %if only 1 player (this should never happen)or not a clean case score just
    %ends
    if(~clean)
        scores_gap(j)=distance_anisotropic(esequence{t}.finalpoints(i,:)',esequence{candidates_t(j)}.finalpoints(candidates(j),:)',trackingparameters.anisotropyvector);
       
        matchings(j,:)=[1,0,0];
        playersstart(j,:)=[i,0,0];
        playersend(j,:)=[candidates(j),0,0];
    end
    %if clean 2 players
    %score 2 options
    if(clean&&isempty(endforwardtrace))
        score1=distance_anisotropic(esequence{t}.finalpoints(startconflictplayers(1),:)',esequence{candidates_t(j)}.finalpoints(endconflictplayers(1),:)',trackingparameters.anisotropyvector);
        score1=score1+distance_anisotropic(esequence{t}.finalpoints(startconflictplayers(2),:)',esequence{candidates_t(j)}.finalpoints(endconflictplayers(2),:)',trackingparameters.anisotropyvector);
        score1=score1/2;
        score2=distance_anisotropic(esequence{t}.finalpoints(startconflictplayers(1),:)',esequence{candidates_t(j)}.finalpoints(endconflictplayers(2),:)',trackingparameters.anisotropyvector);
        score2=score2+distance_anisotropic(esequence{t}.finalpoints(startconflictplayers(2),:)',esequence{candidates_t(j)}.finalpoints(endconflictplayers(1),:)',trackingparameters.anisotropyvector);
        score2=score2/2;
        playersstart(j,:)=[startconflictplayers',0];
        playersend(j,:)=[endconflictplayers',0];
        if(score1<score2)
            scores_gap(j)=score1;
            matchings(j,:)=[1,2,0];
        else
            scores_gap(j)=score2;
            matchings(j,:)=[2,1,0];
        end
        
    end
    
    %if clean 3 players
    %score 3 options
    if(clean&&~isempty(endforwardtrace))
        astart=[startconflictplayers;startbacktrace];
        aend=[endconflictplayers;endforwardtrace];
        playersstart(j,:)=astart;
        playersend(j,:)=aend;
        
        ending_time=candidates_t(j);
        points1=[esequence{t}.finalpoints(astart(1),:);esequence{t}.finalpoints(astart(2),:);esequence{t}.finalpoints(astart(3),:)];
        points2=[esequence{ending_time}.finalpoints(aend(1),:);esequence{ending_time}.finalpoints(aend(2),:);esequence{ending_time}.finalpoints(aend(3),:)];
        [min,config]=scoreTriplePosition(points1,points2,trackingparameters.anisotropyvector(3));
        scores_gap(j)=min/3;
        matchings(j,:)=config;
    end

end
 scores_gap=scores_gap./averagediam;
end

