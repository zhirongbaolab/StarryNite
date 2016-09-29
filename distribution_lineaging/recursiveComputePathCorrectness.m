function esequence=recursiveComputePathCorrectness(esequence,t,i,con,trackingparameters)
%if(isfield(esequence{t},'linkconfidences'))
esequence{t}.path_correctness(i)=con;
if(esequence{t}.suc(i,1)~=-1&esequence{t}.suc_time(i,1)<trackingparameters.endtime)
    [match,matchi,matchj] = checkIfRealSuccessors( t,i,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),esequence );
    esequence=recursiveComputePathCorrectness(...
        esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),con&match,trackingparameters);
    if(esequence{t}.suc(i,2)~=-1&esequence{t}.suc_time(i,2)<trackingparameters.endtime)
        [match,matchi,matchj] = checkIfRealSuccessors( t,i,esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),esequence );
        esequence=recursiveComputePathCorrectness(...
            esequence,esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),con&match,trackingparameters);
    end
end
%end

%{


%pull out wrongly classified cases
confidence_nondiv=posterior(trackingparameters.linkconfidencemodel.nondiv,...
    alllinkconfidencedata.nondiv(:,trackingparameters.linkconfidencemodel.nondivkeep));
cases=alllinkconfidencedata.nondivtrue;

bhc=(~logical(cases)&confidence_nondiv(:,1)<.5);

wronginfo=[];
for i=1:length(bhc)
if(bhc(i))
wronginfo=[wronginfo;linkconfidencedata.nondiv_ti(i,1),...
    esequence{linkconfidencedata.nondiv_ti(i,1)}.finalpoints(linkconfidencedata.nondiv_ti(i,2),:)];
end
end



trackingparameters.endtime=280;
 esequence_con= scoreLinkConfidence(esequence,trackingparameters);

%initialize
for t=1:trackingparameters.endtime
    esequence{t}.path_confidence=zeros(size(esequence{t}.delete));
end
%assign path confidence
startt=50;
for i=1:size(esequence_con{startt}.finalpoints,1)
esequence_con=recursiveComputePathConfidence_min(esequence_con,startt,i,1);
end
esequence_con2=esequence_con;
%delete below .5
for t=50:trackingparameters.endtime-1
esequence_con2{t}.delete(esequence_con2{t}.path_confidence<.9)=true;
end
for t=50:trackingparameters.endtime-1
esequence_con2{t}.suc(logical(esequence_con2{t}.delete),1:2)=-1;
esequence_con2{t}.suc_time(logical(esequence_con2{t}.delete),1:2)=-1;
esequence_con2{t}.pred_time(logical(esequence_con2{t}.delete))=-1;
esequence_con2{t}.pred(logical(esequence_con2{t}.delete))=-1;
end
saveGreedyNucleiFiles(esequence_con2,endtime,outputdirectory,anisotropy,ROIxmin,ROIymin)



startt=20;
for i=1:size(esequence_con{startt}.finalpoints,1)
esequence_con=recursiveComputePathConfidence(esequence_con,startt,i,1);
end
for i=1:size(esequence_con{startt}.finalpoints,1)
esequence_con=recursiveComputePathCorrectness(esequence_con,startt,i,true,trackingparameters);
end

figure
hold on
for t=startt:trackingparameters.endtime-1
good=logical(esequence_con{t}.path_correctness);
    scatter(esequence_con{t}.path_confidence(good)+.05*rand(size(esequence_con{t}.path_confidence(good))),t*ones(size(esequence_con{t}.path_confidence(good))),'b');
  scatter(esequence_con{t}.path_confidence(~good)+.05*rand(size(esequence_con{t}.path_confidence(~good))),t*ones(size(esequence_con{t}.path_confidence(~good))),'r');

end


figure
hold on
for t=startt:trackingparameters.endtime-1
    scatter(esequence_con{t}.path_confidence,t*ones(size(esequence_con{t}.path_confidence)));
end

figure
hold on
for t=1:trackingparameters.endtime-1
    good=logical(esequence_con{t}.;
    scatter(esequence_con{t}.linkconfidences(good),t*ones(size(esequence_con{t}.linkconfidences(good))),'b');
  scatter(esequence_con{t}.linkconfidences(~good),t*ones(size(esequence_con{t}.linkconfidences(~good))),'r');
end



%}