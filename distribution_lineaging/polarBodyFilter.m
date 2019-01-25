function [ esequence] = polarBodyFilter( esequence,starttime,endtime,sizethresh,brightthresh,brightthreshhigh )
%filters cells with diam <sizethresh and mean max slice brightness>brightthresh
for i=starttime:endtime
    for j=1:length(esequence{i}.finaldiams)
        if(esequence{i}.finaldiams(j)<sizethresh&&...
                  esequence{i}.mdiskMax(j)>brightthresh)
%                esequence{i}.maxslice(j)/ esequence{i}.finaldiams(j)>brightthresh)
            esequence{i}.delete(j)=1;
        end
        
    end
end

%second pass if delete is isolated and fails stricter bright threshold
%then undelete it
%otherwise if is connected to something deleted or very bright remains
%deleted
for i=starttime:endtime
    for j=1:length(esequence{i}.finaldiams)
        if(esequence{i}.delete(j))
            deletecurrent=false;
            forwarddelete=i~=length(esequence)&&...
                esequence{i}.suc(j,1)~=-1&&...
                esequence{i+1}.delete(esequence{i}.suc(j,1))==1;
            backdelete=i~=1&&...
                esequence{i}.pred(j)~=-1&&...
                esequence{i-1}.delete(esequence{i}.pred(j))==1; 
            if(esequence{i}.mdiskMax(j)<brightthreshhigh&&...
                    ~forwarddelete&&~backdelete)
                esequence{i}.delete(j)=0;
            else
                %if deleted make sure its not linking itself to anything
                if(esequence{i}.suc(j,1)~=-1)
                    %unlink if from suc
                    esequence{i+1}.pred(esequence{i}.suc(j,1))=-1;
                    esequence{i}.suc(j,1)=-1;
                end
                 if(esequence{i}.pred(j)~=-1)
                    %unlink if from pred
                    esequence{i-1}.suc(esequence{i}.pred(j),1)=-1;
                    esequence{i}.pred(j)=-1;
                end
            end
        end
        
    end
end


end

