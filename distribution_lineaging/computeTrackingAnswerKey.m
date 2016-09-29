function [esequence]=computeTrackingAnswerKey(esequence,endtime)
%computeTrackingAnswerKey computes correct tracking correspndences for detected points
%given matches with corrected answer and translated (from cell tag numbers to row indicies)
%corrected answer matches i.e. computes exactly what the tracking alg
%should give you for each detected nucleus
%-1 if FP
%nuc index of next frame sucessor if present
%nuc index  and time index of next detected nucleus in sequence if FN
%occurs ahead

for time=1:endtime


    if (time<endtime) %last timepoint is all -1
        %initialize answer        
        esequence{time}.correct_suc=zeros(size(esequence{time}.finalpoints,1),2);
        esequence{time}.correct_suc_time=zeros(size(esequence{time}.finalpoints,1),2);
        esequence{time}.correct_pred=-1*ones(size(esequence{time}.finalpoints,1),1);
        esequence{time}.correct_pred_time=-1*ones(size(esequence{time}.finalpoints,1),1);
   
    %iterate over detected points
   e=esequence{time};
   for j=1:size(e.finalpoints,1)
       %this point is matched (not FP)|| not  a death
       for child=1:2
           if (~(e.matches(j)==-1)&&~(e.corrected_sucessors(e.matches(j),child)==-1))
               sucessor_match=esequence{time+1}.matchesr(e.corrected_sucessors(e.matches(j),child));
               if(sucessor_match~=-1) % the corrected sucessor to match is also matched
                   esequence{time}.correct_suc(j,child)=sucessor_match;
                   esequence{time}.correct_suc_time(j,child)=time+1;
               else %sucessor is FN look ahead for next sucessor that exists link to that
                   unmatched=true;
                   dead=false;
                   match=-1;
                   matchtime=-1;
                   localtime=time;
                   localsuccessor=e.corrected_sucessors(e.matches(j),child);
                   dead=localtime>=endtime-1;
                   while(unmatched&&~dead)
                       %it is not matched search ahead for match till find it or
                       %hit a death
                       if (localsuccessor==-1)
                           dead=true;
                       else
                           localsucessor_match=esequence{localtime+1}.matchesr(localsuccessor);
                           
                           if(localsucessor_match~=-1)
                               match=localsucessor_match;
                               matchtime=localtime+1;
                               unmatched=false;
                           end
                           localsuccessor=esequence{localtime+1}.corrected_sucessors(localsuccessor,1);                  
                           localtime=localtime+1;
                           if (localtime==endtime-1)
                               dead=true;
                           end
                       end
                   end
                   esequence{time}.correct_suc(j,child)=match;
                   esequence{time}.correct_suc_time(j,child)=matchtime;
               end
           else %FP has no match || is a death
               esequence{time}.correct_suc(j,child)=-1;
               esequence{time}.correct_suc_time(j,child)=-1;
           end
       end
   end
    else
          %initialize answer  
           esequence{time}.correct_suc=-1*ones(size(esequence{time}.finalpoints,1),2);
    esequence{time}.correct_suc_time=-1*ones(size(esequence{time}.finalpoints,1),2);

    end
end
%copy correct suc into correct pred
for time=1:endtime-1
    
    %iterate over detected points
    e=esequence{time};
    for j=1:size(e.finalpoints,1)
        if(e.correct_suc(j,1)~=-1)
        esequence{e.correct_suc_time(j,1)}.correct_pred(e.correct_suc(j,1))=j;
        esequence{e.correct_suc_time(j,1)}.correct_pred_time(e.correct_suc(j,1))=time;
        end
        if(e.correct_suc(j,2)~=-1)
        esequence{e.correct_suc_time(j,2)}.correct_pred(e.correct_suc(j,2))=j;
        esequence{e.correct_suc_time(j,2)}.correct_pred_time(e.correct_suc(j,2))=time;
        end
    end
  
end

