function [dist ] = nnInterfeerenceConfidenceFeature( esequence,t,i,direction,trackingparameters )
%compute at time t pos i the vicinity of a
%end if direction=0
%start if direction =1
%normalize this by local nn distance at i within timepoint
justone=false;
if(t>0&&t<length(esequence))
    suc=esequence{t}.suc(i,:);
    suct=esequence{t}.suc_time(i,:);
    
    if (direction==0)
        interesting=esequence{t}.sucessor_suitors{i};
        justone=length(esequence{t})==1;
        interesting=interesting(interesting~=esequence{t}.suc(i,1)&interesting~=esequence{t}.suc(i,2));
    else
        justone=length(esequence{suct(1)})==1;
        interesting=esequence{suct(1)}.predecessor_suitors{suc(1)};
        if(suc(2)~=-1)
            interesting=[interesting;esequence{suct(2)}.predecessor_suitors{suc(2)}];
        end
        interesting=unique(interesting);
        interesting=interesting(interesting~=i);
        
    end
else
    interesting=[];
end

if(isempty(interesting)||justone)
    dist=[10,10];%10xnn sufficient%  dist=max(esequence{t}.selfdistance);%if there are no ends use dummy value of most distant point
else
    linkdistxy=distance(esequence{t}.finalpoints(i,1:2)',esequence{suct(1)}.finalpoints(suc(1),1:2)');
    linkdistz=distance(esequence{t}.finalpoints(i,3),esequence{suct(1)}.finalpoints(suc(1),3))*trackingparameters.anisotropyvector(3);
    if(suc(2)~=-1)
        linkdistxy=max(linkdistxy,distance(esequence{t}.finalpoints(i,1:2)',esequence{suct(2)}.finalpoints(suc(2),1:2)'));
        linkdistz=max(linkdistz,distance(esequence{t}.finalpoints(i,3),esequence{suct(2)}.finalpoints(suc(2),3))*trackingparameters.anisotropyvector(3));
    end
    if (direction==0)
        %note that this depends on divisions never being to different times
        %for correctness
        [mdist,mini]=min(distance_anisotropic(esequence{t}.finalpoints(i,:)',...
            esequence{suct(1)}.finalpoints(interesting,:)',trackingparameters.anisotropyvector));
        minxy=distance(esequence{t}.finalpoints(i,1:2)',esequence{suct(1)}.finalpoints(mini,1:2)');
        minz=distance(esequence{t}.finalpoints(i,3),esequence{suct(1)}.finalpoints(mini,3))*trackingparameters.anisotropyvector(3);
    else
        [mdist,mini]=min(distance_anisotropic(esequence{suct(1)}.finalpoints(suc(1),:)',...
            esequence{t}.finalpoints(interesting,:)',trackingparameters.anisotropyvector));
        s2=false;
        if(suc(2)~=-1)
            [mdist2,mini2]=min(distance_anisotropic(esequence{suct(2)}.finalpoints(suc(2),:)',...
                esequence{t}.finalpoints(interesting,:)',trackingparameters.anisotropyvector));
            
            if (mini2<mini)
                mdist=mdist2;
                s2=true;
            end
        end
        if(s2)
            minxy=distance(esequence{t}.finalpoints(mini,1:2)',esequence{suct(2)}.finalpoints(suc(2),1:2)');
            minz=distance(esequence{t}.finalpoints(mini,3),esequence{suct(2)}.finalpoints(suc(2),3))*trackingparameters.anisotropyvector(3);
        else
            minxy=distance(esequence{t}.finalpoints(mini,1:2)',esequence{suct(1)}.finalpoints(suc(1),1:2)');
            minz=distance(esequence{t}.finalpoints(mini,3),esequence{suct(1)}.finalpoints(suc(1),3))*trackingparameters.anisotropyvector(3);
            
        end
    end
    dist=[linkdistxy/minxy,linkdistz/minz];
end
  
end

