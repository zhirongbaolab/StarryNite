function [ data ] = calculateWideWindowCellPairVectorNondivision(esequence,t,i,t2,j,anisotropyvector,wideWindow)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gapvector=calculateCellPairVectorNondivision_wdiam(esequence{t},i,esequence{t2},j,anisotropyvector);
%calculateWideWindowCellPairVectorNondivision(esequence,t,i,esequence{t}.suc_time(i,d),esequence{t}.suc(i,d),trackingparameters.anisotropyvector,trackingparameters.wideWindow)
beforevectors=[];
curr=i;
curr_t=t;
pred=esequence{t}.pred(i);
pred_t=esequence{t}.pred_time(i);
count=1;
%while track doesnt end,within window and no division 
while (pred~=-1&&count<=wideWindow&&esequence{pred_t}.suc(pred,2)==-1)
    beforevectors=[beforevectors;...
        calculateCellPairVectorNondivision_wdiam(...
        esequence{pred_t},pred,esequence{curr_t},curr,anisotropyvector)...
        ];
    curr=pred;
    curr_t=pred_t;
    bkpred=pred;
    pred=esequence{pred_t}.pred(pred);
    pred_t=esequence{pred_t}.pred_time(bkpred);
    count=count+1;    
end
aftervectors=[];
curr=j;
curr_t=t2;
suc=esequence{t2}.suc(j,1);
suc_t=esequence{t2}.suc_time(j,1);
count=1;
%while track doesnt end,within window and no division 
while (suc~=-1&&count<=wideWindow&&esequence{curr_t}.suc(curr,2)==-1)
    aftervectors=[aftervectors;...
        calculateCellPairVectorNondivision_wdiam(...
        esequence{curr_t},curr,esequence{suc_t},suc,anisotropyvector)...
        ];
    curr=suc;
    curr_t=suc_t;
    suc=esequence{suc_t}.suc(curr,1);
    suc_t=esequence{suc_t}.suc_time(curr,1);
    count=count+1;    
end
vecs=[beforevectors;aftervectors];
if(size(vecs,1)>1)
    vecs=mean(vecs);
end
if(isempty(vecs))
    vecs=1*ones(1,5);
end

data=gapvector./vecs;
data(isnan(data)|isinf(data))=0;
%if(~isempty(find(isnan(data))))
%    'odd'
%end
end

