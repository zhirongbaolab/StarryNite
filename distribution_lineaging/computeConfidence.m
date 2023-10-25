function [convector,consum]=computeConfidence(esequence,t,i,minsize)
%integrate confidence and and confidence vector over a specified number of
%frames
convector=[];
consum=0;
for j=1:minsize
convector=[convector;esequence{t}.confidencevector(i,:)];
consum=consum+esequence{t}.confidences(i,:);
ttemp=t;
t=esequence{t}.suc_time(i,1);
i=esequence{ttemp}.suc(i,1);
    if i==-1
        break
    end
end

end

