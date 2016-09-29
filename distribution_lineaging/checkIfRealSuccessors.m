function [match,matchi,matchj ] = checkIfRealSuccessors( ti,is,tj,js,esequence )
% given a list of nuclearindicies at time ti and indicies at tj
%computes whether some j is a sucessor to some i according to the
%cannonical answer
match=false;
matchi=[];
matchj=[];
for j=1:length(js)
    for i=1:length(is)
        current=js(j);
        currentt=tj;
        dead=false;
        while (currentt>ti&~match&~dead)
            tempcur=current;
            current=esequence{currentt}.correct_pred(current);
            currentt=esequence{currentt}.correct_pred_time(tempcur);
            if(current==-1)
                dead=true;
            end
            
            if(currentt==ti)
                if(current==is(i))
                    match=true;
                    matchi=i;
                    matchj=j;
                end
            end
            
        end
    end
end

end

