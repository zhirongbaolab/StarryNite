function [convector,consum]=computeConfidenceback(esequence,t,i,minsize)
%integrate confidence and and confidence vector over a specified number of
%frames
convector=[];
consum=0;
for j=1:minsize
convector=[convector;esequence{t}.confidencevector(i,:)];
consum=consum+esequence{t}.confidences(i,:);
ttemp=t;
t=esequence{t}.pred_time(i);
i=esequence{ttemp}.pred(i);
end

end

