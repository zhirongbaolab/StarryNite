function good=filterCandidateLength(suitors,deletemarked)
%iterate over suitors and return length of those which are not labled as FP
good=0;
for i=1:length(suitors)
    if ~deletemarked(suitors(i))
        good=good+1;
    end
end

end

