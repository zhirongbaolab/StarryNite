function [error ,correctflags,computedclasses,model] = ...
    compute_classificationerror( classtags, alldaughterdata,allbackdata,...
    allforwarddata,trulyambigious,FullyDivLooking,...
    FNDivLooking,FullyFPLooking,DirtyFPLooking,DivFPLooking)
%train a classifier model based on data and labels and return the model and its errors
distributiontypemain='kernel';
correct=0;
correctflags=ones(size(classtags));
computedclasses=zeros(size(classtags));
model=[];

%ambigious

data=[alldaughterdata(trulyambigious,:),allbackdata(trulyambigious,:),allforwarddata(trulyambigious,:)];

if(length(find(trulyambigious))>1)
    
classes=classtags(trulyambigious);
%'hack'
%classes(classes==3)=2;
%data(isinf(data))=10;
testclass=NaiveBayes.fit(data,classes,'distribution',distributiontypemain);
%testclass=classregtree(data,classes,'Method','classification');
testpred=predict(testclass,data);
%testpred=str2num(cell2mat(testclass.eval(data,10)));
model.ambigious=testclass;

testconfusion=confusionmat(classes,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=testpred~=classes;

casetargets=find(trulyambigious);
correctflags(casetargets(wrong))=0;
computedclasses(casetargets)=testpred;

else
    'null model for ambigious'
    %in case of no data make null model with random data and other label
    %for all
    model.ambigious=NaiveBayes.fit(rand(20,size(data,2)),zeros(1,20),'distribution',distributiontypemain);
end

%fp or div (both just bif data)

data=[alldaughterdata(FullyFPLooking|FullyDivLooking,:)];

if(length(find(FullyFPLooking|FullyDivLooking))>1)
classes=classtags(FullyFPLooking|FullyDivLooking);
classes(classes==2)=0; %fn is invalid choice
%data(isinf(data))=10;
testclass=NaiveBayes.fit(data,classes,'distribution',distributiontypemain);
%testclass=classregtree(data,classes,'Method','classification');
testpred=predict(testclass,data);
%testpred=str2num(cell2mat(testclass.eval(data,10)));

model.fp_div=testclass;

testconfusion=confusionmat(classes,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=testpred~=classes;
casetargets=find(FullyFPLooking|FullyDivLooking);
correctflags(casetargets(wrong))=0;
computedclasses(casetargets)=testpred;
else
    'null model for fp_div'
    %in case of no data make null model with random data and other label
    %for all
    model.fp_div=NaiveBayes.fit(rand(20,size(data,2)),zeros(1,20),'distribution',distributiontypemain);
end

%dirty Fp or FN div
data=[alldaughterdata(DirtyFPLooking|FNDivLooking,:),allbackdata(DirtyFPLooking|FNDivLooking,:)];

if(length(find(DirtyFPLooking|FNDivLooking))>1)
classes=classtags(DirtyFPLooking|FNDivLooking);
%data=data(classes~=0,:);
%'hack'
%classes(classes==3)=0;

%classes=classes(classes~=0);
%all are valid but remove singleton other class
%data(isinf(data))=10;
testclass=NaiveBayes.fit(data,classes,'distribution',distributiontypemain);
%testclass=classregtree(data,classes,'Method','classification');
testpred=predict(testclass,data);
%testpred=str2num(cell2mat(testclass.eval(data,10)));
model.dirtyfp_fn=testclass;

testconfusion=confusionmat(classes,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=testpred~=classes;
casetargets=find(DirtyFPLooking|FNDivLooking);
%casetargets=find(DirtyFPLooking);
correctflags(casetargets(wrong))=0;
computedclasses(casetargets)=testpred;
else
    'null model for dirtyfp'
    %in case of no data make null model with random data and other label
    %for all
    model.dirtyfp_fn=NaiveBayes.fit(rand(20,size(data,2)),zeros(1,20),'distribution',distributiontypemain);
end

%div Fp %div
data=[alldaughterdata(DivFPLooking,:),allforwarddata(DivFPLooking,:)];

if(length(find(DivFPLooking))>1)
classes=classtags(DivFPLooking);
classes(classes~=1&classes~=3)=0; %only div, FP and other are valid outcomes

%data(isinf(data))=10;
testclass=NaiveBayes.fit(data,classes,'distribution',distributiontypemain);
%testclass=classregtree(data,classes,'Method','classification');
testpred=predict(testclass,data);
%testpred=str2num(cell2mat(testclass.eval(data,10)));
model.divfp=testclass;

testconfusion=confusionmat(classes,testpred)
for i=1:size(testconfusion,1)
    correct=correct+testconfusion(i,i);
end
wrong=find(testpred~=classes);
casetargets=find(DivFPLooking);
correctflags(casetargets(wrong))=0;
computedclasses(casetargets)=testpred;

else
    'null model for divfp'
    %in case of no data make null model with random data and other label
    %for all
    model.divfp=NaiveBayes.fit(rand(20,size(data,2)),zeros(1,20),'distribution',distributiontypemain);
end

error=size(alldaughterdata,1)-correct;

end


