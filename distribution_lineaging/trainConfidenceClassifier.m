function [model]=trainConfidenceClassifier(confidencedata)
model=[];
correct=0;
modeltype='normal';
dimreduction=false;

keepers=ones(1,size(confidencedata.div,2));

%in normal mode drop all zero variance features
if(strcmp(modeltype,'normal'))
    sigma=std(confidencedata.div(logical(confidencedata.divtrue),:));
    gsize=length(find(confidencedata.divtrue));
    badCols = sigma <= gsize * eps(max(sigma));
      sigma=std(confidencedata.div(~logical(confidencedata.divtrue),:));
    gsize=length(find(~confidencedata.divtrue));
    badCols = badCols|(sigma <= gsize * eps(max(sigma)));
    keepers(badCols)=false;
end
if (dimreduction)
for k=1:length(keepers)
%div
%with feature k
testclassdiv=NaiveBayes.fit(confidencedata.div(:,logical(keepers)),confidencedata.divtrue...
    ,'distribution',modeltype);
testpred=predict(testclassdiv,confidencedata.div(:,logical(keepers)));
wrong=length(find(testpred~=confidencedata.divtrue));
%without feature k
keepers(k)=0;
testclassdiv=NaiveBayes.fit(confidencedata.div(:,logical(keepers)),confidencedata.divtrue...
    ,'distribution',modeltype);
testpred=predict(testclassdiv,confidencedata.div(:,logical(keepers)));
wrong2=length(find(testpred~=confidencedata.divtrue));
if(wrong2>wrong)
    keepers(k)=1;
end
end
end
%Div final plass with ultimate model
testclassdiv=NaiveBayes.fit(confidencedata.div(:,logical(keepers)),confidencedata.divtrue...
    ,'distribution',modeltype);
testpred=predict(testclassdiv,confidencedata.div(:,logical(keepers)));
testconfusion=confusionmat(confidencedata.divtrue,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=testpred~=confidencedata.divtrue;
model.div=testclassdiv;
model.divkeep=logical(keepers);



%gap train
keepers=ones(1,size(confidencedata.gap,2));
keepers=logical(keepers);
%in normal mode drop all zero variance features
if(strcmp(modeltype,'normal'))
    sigma=std(confidencedata.gap(logical(confidencedata.gaptrue),:));
    gsize=length(find(confidencedata.gaptrue));
    badCols = sigma <= gsize * eps(max(sigma));
      sigma=std(confidencedata.gap(~logical(confidencedata.gaptrue),:));
    gsize=length(find(~confidencedata.gaptrue));
    badCols = badCols|(sigma <= gsize * eps(max(sigma)));
    keepers(badCols)=false;
end
if(dimreduction)
for k=1:length(keepers)
testclassgap=NaiveBayes.fit(confidencedata.gap(:,keepers),confidencedata.gaptrue...
    ,'distribution',modeltype);
testpred=predict(testclassgap,confidencedata.gap(:,keepers));
wrong=length(find(testpred~=confidencedata.gaptrue));
keepers(k)=false;
testclassgap=NaiveBayes.fit(confidencedata.gap(:,keepers),confidencedata.gaptrue...
    ,'distribution',modeltype);
testpred=predict(testclassgap,confidencedata.gap(:,keepers));
wrong2=length(find(testpred~=confidencedata.gaptrue));
if(wrong2>wrong)
    keepers(k)=1;
end
end
end
%gap final test pass
testclassgap=NaiveBayes.fit(confidencedata.gap(:,keepers),confidencedata.gaptrue...
    ,'distribution',modeltype);
testpred=predict(testclassgap,confidencedata.gap(:,keepers));
confidence=posterior(testclassgap,confidencedata.gap(:,keepers));
testconfusion=confusionmat(confidencedata.gaptrue,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=length(find(testpred~=confidencedata.gaptrue));
model.gap=testclassgap;
model.gapkeep=logical(keepers);
%{
bins=linspace(0,1,10);
    goodc=histc(confidence(logical(confidencedata.gaptrue),1),bins);
      badc=histc(confidence(~logical(confidencedata.gaptrue),1),bins);
%}



%train nondiv
keepers=ones(1,size(confidencedata.nondiv,2));
keepers=logical(keepers);
%in normal mode drop all zero variance features
if(strcmp(modeltype,'normal'))
    sigma=std(confidencedata.nondiv(logical(confidencedata.nondivtrue),:));
    gsize=length(find(confidencedata.nondivtrue));
    badCols = sigma <= gsize * eps(max(sigma));
      sigma=std(confidencedata.nondiv(~logical(confidencedata.nondivtrue),:));
    gsize=length(find(~confidencedata.nondivtrue));
    badCols = badCols|(sigma <= gsize * eps(max(sigma)));
    keepers(badCols)=false;
end
%'skipping keep training on non division'

if(dimreduction)
for k=1:length(keepers)

testclassnondiv=NaiveBayes.fit(confidencedata.nondiv(:,keepers),confidencedata.nondivtrue...
    ,'distribution',modeltype);
testpred=predict(testclassnondiv,confidencedata.nondiv(:,keepers));
wrong=length(find(testpred~=confidencedata.nondivtrue));
keepers(k)=false;
testclassnondiv=NaiveBayes.fit(confidencedata.nondiv(:,keepers),confidencedata.nondivtrue...
    ,'distribution',modeltype);
testpred=predict(testclassnondiv,confidencedata.nondiv(:,keepers));
wrong2=length(find(testpred~=confidencedata.nondivtrue));

if(wrong2>wrong)
    keepers(k)=true;
end
end
end
%}
%final test pass nondiv
testclassnondiv=NaiveBayes.fit(confidencedata.nondiv(:,keepers),confidencedata.nondivtrue...
    ,'distribution',modeltype);
testpred=predict(testclassnondiv,confidencedata.nondiv(:,keepers));
testconfusion=confusionmat(confidencedata.nondivtrue,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=testpred~=confidencedata.nondivtrue;
model.nondiv=testclassnondiv;
model.nondivkeep=logical(keepers);


%create a false real death to keep training from crashing though trained
%model will be useless
if(length(find(confidencedata.deathtrue))<2)
    'creating dummy death because no real deaths found'
    confidencedata.deathtrue(1:2)=true;
end

keepers=ones(1,size(confidencedata.death,2));
keepers=logical(keepers);

%in normal mode drop all zero variance features
if(strcmp(modeltype,'normal'))
    sigma=std(confidencedata.death(logical(confidencedata.deathtrue),:));
    gsize=length(find(confidencedata.deathtrue));
    badCols = sigma <= gsize * eps(max(sigma));
      sigma=std(confidencedata.death(~logical(confidencedata.deathtrue),:));
    gsize=length(find(~confidencedata.deathtrue));
    badCols = badCols|(sigma <= gsize * eps(max(sigma)));
    keepers(badCols)=false;
end
if(dimreduction)
for k=1:length(keepers)  
    testclassdeath=NaiveBayes.fit(confidencedata.death(:,keepers),confidencedata.deathtrue...
        ,'distribution',modeltype);
    testpred=predict(testclassdeath,confidencedata.death(:,keepers));
    wrong=length(find(testpred~=confidencedata.deathtrue));
    testclassdeath=NaiveBayes.fit(confidencedata.death(:,keepers),confidencedata.deathtrue...
        ,'distribution',modeltype);
    testpred=predict(testclassdeath,confidencedata.death(:,keepers));
    wrong2=length(find(testpred~=confidencedata.deathtrue));
    
    if(wrong2>wrong)
        keepers(k)=true;
    end
end
end
testclassdeath=NaiveBayes.fit(confidencedata.death(:,keepers),confidencedata.deathtrue...
    ,'distribution',modeltype);
testpred=predict(testclassdeath,confidencedata.death(:,keepers));
testconfusion=confusionmat(confidencedata.deathtrue,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=testpred~=confidencedata.deathtrue;
model.death=testclassdeath;
model.deathkeep=logical(keepers);
%}

