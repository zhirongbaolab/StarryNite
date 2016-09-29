function [sscore,mscore,range1,range2]=mergeDecision(n1,n2,n1c,n2c,centers,ranges,logodds,xymax,xymaximavals,xydetdiameters,anisotropy,xycoverage)

%for moment this is degenerating into is mergescore >0
mscore=mergeScore2(n1,n2,n1c,n2c,centers,ranges,xymax,xymaximavals,xydetdiameters,anisotropy,xycoverage);

[sscore,range1,range2]=splitScore(n1,n2,n1c,n2c,centers,ranges,logodds);
%now returning split ranges and 2 scores regardless
%{
if(mscore>sscore)
    %merge=true;
     range1=ranges{n1};
   range2=ranges{n2};
else
   %use the ranges that come from splitscore for each if split
   
end
%}