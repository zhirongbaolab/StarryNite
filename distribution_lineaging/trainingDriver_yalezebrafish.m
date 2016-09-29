%{
%remember you cano nly run this for the embryo that is currently in
% workspace and need to reload it for others to get points, dumbass

emb=lin;
cases=allbifurcationinfo(emb).refclassificationvector==0;
cases=allbifurcationinfo(emb).refclassificationvector~=allbifurcationinfo(emb).computedclassificationvector;
cases=allbifurcationinfo(emb).refclassificationvector==0&allbifurcationinfo(emb).classround'==2;

info=[];
for i=1:length(cases)
    if(cases(i))
        info=[info;allbifurcationinfo(emb).removed(i,1),...
            esequence{allbifurcationinfo(emb).removed(i,1)}.finalpoints(...
        allbifurcationinfo(emb).removed(i,2),:)];
    end
end
if(ROI)
info(:,2)=info(:,2)+ROIxmin;
info(:,3)=info(:,3)+ROIymin
end
  info=[info,allbifurcationinfo(emb).refclassificationvector(cases)',allbifurcationinfo(emb).computedclassificationvector(cases)'];

cmats={};
for i=1:length(allbifurcationinfo)
cmats{i}=confusionmat(allbifurcationinfo(i).refclassificationvector,allbifurcationinfo(i).computedclassificationvector);
end
%}
%evaluate error over set of lineages for marker 1

basematdir='L:\santella\Hobert_sox_data\20120427_otls333_dualcolor\';
basedir='L:\santella\Hobert_sox_data\20120427_otls333_dualcolor\20120427_otls333_dualcolor\';
lineages={... 
'combined_emb_unflipgood1',...
};
lineagedir={'20120427_otls333_dualcolor'};
lineageimage={'combined'};


basematdir='L:\santella\yale_zebrafish\';
basedir='L:\santella\yale_zebrafish\annot\yale_zebrafish\';
lineages={'yale_zebrafish_emb_training_subset1'};
matlineages={'yale_zebrafish_emb_baseline5min_nobif5_8p75rangethresh_mergetrain1'};
%first version was both lineages
%lineages={'yale_zebrafish_emb_newparam_retrainmodel_kernel_21minwholetest1'};
%matlineages={'yale_zebrafish_emb_newparam_retrainmodel_kernel_21minwholetest1'};

lineagedir={'yale_zebrafish'}
lineageimage={'yale_zebrafish'};
eightcelltimes=30*ones(size(lineagedir));
%eightcelltimes=[35,35,25,25,35,25,35,35,35,35];
eightcelltimes=eightcelltimes+20;%go to AB8 instead of P08
nondivthress=[];

edittimes=[10];

%red 5 training better for doing confidence training because lack 'wrong'
%polar links
%when these are used for training filtering should be off
%{
lineages={'ZD_BV82_WT_20100809_2_s1_emb_linbiftest_lowthresh1','ZD_RW10425_WT_20100412_2_s1_emb1',...
    'ZD_RW10425_WT_20100412_2_s1_emb2','ZD_RW10425_WT_20100412_2_s1_emb3','ZD_RW10434_WT_20110429_2_s1_emb_linbiftest1'};
lineagedir={'ZD_BV82_WT_20100809_2','ZD_RW10425_WT_20100412_2',...
    'ZD_RW10425_WT_20100412_2','ZD_RW10425_WT_20100412_2','ZD_RW10434_WT_20110429_2'};
lineageimage={'ZD_BV82_WT_20100809_2_s1','ZD_RW10425_WT_20100412_2_s1',...
    'ZD_RW10425_WT_20100412_2_s1','ZD_RW10425_WT_20100412_2_s1','ZD_RW10434_WT_20110429_2_s1'};
edittimes=[200,192,201,185,185]
eightcelltimes=[20,20,30,20,20];
basedir='L:\santella\unzipped_lineages\test_data\'
%}





%'green mode'
%parameterConfigurationGreen
'red mode'
parameterConfiguration

%trainingmode refers to bif. classifier training not confidence
%uses oracle
trackingparameters.trainingmode=true;
trackingparameters.recordanswers=true;
answerkey=true;
outputtrimmed=false;
errors=cell(1,length(lineages));
'note set to not replace division model because no divisions to create it from'
replacemodel=false;
allexpectedchange={};
evalforced=false;
evalfinal=false;
allbifurcationinfo=[];

for lin=1:length(lineages)
    endtime=edittimes(lin);
    trackingparameters.endtime=endtime;
    %'loader for files in zhuos directory'
     %load([basedir,lineagedir{lin},'/',lineages{lin},'_fullmatlabresult.mat']);
    'loader for training files'
     load([basematdir,matlineages{lin},'_fullmatlabresult.mat']);
    trackingparameters.anisotropyvector=[1,1,anisotropy];
    parameters.anisotropyvector=[1,1,anisotropy];
    
    embryonumbers = {};
    nucleidir=basedir;
    embryonumbers_c={[lineages{lin},'_edited\nuclei\']};
    %endtime=edittimes(lin)+11;
    outputdirectory=[lineages{lin},'/nuclei/'];
    if (answerkey)
        train_tracking_statistics_function;
        train_confidence_function;
    end
    endtime=endtime+10;
    trackingparameters.endtime=trackingparameters.endtime+10;
    
    %for live run reset to all timepoints
%   'process all time points for live trim'
%    
%   'process 250 for top down'  
% trackingparameters.endtime=min(250,length(esequence)); 
% trackingparameters.endtime=length(esequence);
 %       endtime=trackingparameters.endtime;
    
  
    
    trackingparameters.recordanswers=false;
    
    tracking_driver_new_classifier_based_version;
    
    %time adjustment back for final evaluation
    %endtime=endtime-9;
    %trackingparameters.endtime=trackingparameters.endtime-9;
    %reset back for final evaluation
    if(trackingparameters.trainingmode==false)
        endtime=endtime-11;
        trackingparameters.endtime=trackingparameters.endtime-11;
    end
    
    
    %section for storing training info
    ncells=[];
    for i=1:size(removed,1)
        % ncells=[ncells,length(esequence{removed(i,1)}.FP)];
        ncells=[ncells,mean(esequence{removed(i,1)}.selfdistance)./mean(esequence{removed(i,1)}.finaldiams)];
    end
    if (trackingparameters.trainingmode) %training mode reall means training for bifurcation model
        allbifurcationinfo(lin).Divdata=Divdata;
        allbifurcationinfo(lin).Tripledata=Tripledata;
        allbifurcationinfo(lin).NoDivdata=NoDivdata;
        allbifurcationinfo(lin).ncells=ncells;
        allbifurcationinfo(lin).removed=removed;
        allbifurcationinfo(lin).confidenceData=confidenceData;
        allbifurcationinfo(lin).splitFNMatchScore=splitFNMatchScore;
        allbifurcationinfo(lin).BifurcationMeasures=BifurcationMeasures;
        allbifurcationinfo(lin).classround=classround;
        allbifurcationinfo(lin).computedclassificationvector=computedclassificationvector;
        allbifurcationinfo(lin).refclassificationvector= refclassificationvector;
    end
    if(~outputtrimmed)
        [ linkconfidencedata] ...
            = extractTrainingConfidenceDataVectors( esequence,trackingparameters,embryonumbers_c,nucleidir,ROI,ROIxmin,ROIymin  );
        allbifurcationinfo(lin).linkconfidencedata=linkconfidencedata;
    end
    if(evalfinal)
        
        mkdir(outputdirectory);
        saveGreedyNucleiFiles(esequence,endtime,outputdirectory,anisotropy,ROIxmin,ROIymin);
        zipname=[lineages{lin},'/',embryonumber,'_',suffix,'.zip'];
        zip(zipname,[outputdirectory,'']);
        
        uneddir=outputdirectory;
        eddir=[basedir,lineages{lin},'_edited\nuclei\'];
        %    uneddir=[basedir,lineages{i},'\nuclei\'];
        test=evaluate_lineage_error(eddir,uneddir,endtime,anisotropy);
        test{37}=test{10}+test{21}+test{32};
        test{38}=test{12}+test{23}++test{34};
        test{39}=sum(test{37}+test{38});
        test{40}=sum(test{37});
        test{41}=sum(test{38});
        
        errors{lin}=test;
    end
    if(trackingparameters.trainingmode)
        allexpectedchange{lin}=expected_corrections;
    end
    
    
    if(outputtrimmed)
        esequence_con= scoreLinkConfidence(esequence,trackingparameters);
        %{ 
        %cutting trimming version
        for t=1:trackingparameters.endtime-1
            for i=1:length(esequence_con{t}.suc)
                if(esequence_con{t}.linkconfidences(i)<.8)
                    suct=esequence_con{t}.suc_time(i,:);
                    suc=esequence_con{t}.suc(i,:);
                     if(suc(1)~=-1)
                    %wipe out pred of suc
                    esequence_con{suct(1)}.pred(suc(1))=-1;
                    esequence_con{suct(1)}.pred_time(suc(1))=-1;
                     end
                    if(suc(2)~=-1)
                        esequence_con{suct(2)}.pred(suc(2))=-1;
                        esequence_con{suct(2)}.pred_time(suc(2))=-1;
                    end
                    esequence_con{t}.suc(i,:)=[-1,-1];
                    esequence_con{t}.suc_time(i,:)=[-1,-1];
                end
            end
        end
        %}
        %clipped trimming version
      %{
        %initialize
        for t=1:trackingparameters.endtime
            esequence_con{t}.path_confidence=zeros(size(esequence{t}.delete));
        end
        %assign path confidence
        startt=eightcelltimes(lin);
        for i=1:size(esequence_con{startt}.finalpoints,1)
            esequence_con=recursiveComputePathConfidence_min(esequence_con,startt,i,1);
        end
        %delete below .8 previously
        for t=startt:trackingparameters.endtime
            esequence_con{t}.delete(esequence_con{t}.path_confidence<.8)=true;
        end
      %}
        %output
         % finaloutputdirectory='l:/santella/lineage_automerge/80version_liberalmodel_ab8_newversiontrim/';
        finaloutputdirectory='L:\santella\lineage_automerge\80version_conservativemodel_ab8_wholescore\';
        trimmed_outputdirectory=['nuclei/'];
        mkdir(trimmed_outputdirectory);
        saveGreedyNucleiFilesAndConfidence(esequence_con,trackingparameters.endtime, trimmed_outputdirectory,anisotropy,ROIxmin,ROIymin)
        zipname=['trim_',embryonumber,'_',suffix,'.zip'];
        
        zip(zipname, trimmed_outputdirectory);
        rmdir( trimmed_outputdirectory,'s');
        
      
                          
        movefile (zipname, finaloutputdirectory)
        outputXMLfile([finaloutputdirectory,lineages{lin},'_edited.xml'],...
            xyres,zres,slices,[finaloutputdirectory,zipname],...
            ['l:/duz/project/imaging/',lineagedir{lin},'/'],...
            lineageimage{lin},false);
        %}
          %save result as mat file
       save([finaloutputdirectory,lineages{lin},'.mat']);
    end
     finaloutputdirectory='L:\santella\lineage_automerge\80version_conservativemodel_ab8_wholescore\';
        save([finaloutputdirectory,lineages{lin},'.mat']);
   
    
end
%{
if(trackingparameters.trainingmode)
    multiple_embryo_train;
else
    multiple_embryo_train_confidence;
end

%}
    
           