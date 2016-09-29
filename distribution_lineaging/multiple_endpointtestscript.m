%script used in testing multiple thresholds for noneasy 1-1, 
%not confirmed good in light of revisions but should be
errors={};
count=1;
evalforced=false;
basedir='L:\santella\unzipped_lineages\';
basedir='L:\santella\unzipped_lineages\training\'

nucleidir=basedir;
replacemodel=true;
%eddir='G:\My Documents\latetest\edited_nuclei\'
%embryonumber='ZD_RW10348_PAL-1_20110318_2_s3_emb1';
embryonumber='ZD_BV82_PIE-1_20110315_1_s2_emb1';
embryonumber='ZD_BV82_APX-1_20110415_1_s2_emb1';
embryonumbers_c={};
embryonumbers_c{1}=[embryonumber,'_edited\nuclei'];
lin=1;
nondivthress=[];

eddir=[basedir,embryonumber,'_edited\nuclei\'];
outputdirectory=[embryonumber,'./nuclei/'];
  load([basedir,embryonumber,'_fullmatlabresult.mat']);
  endtime=190;
   train_tracking_statistics_function;
 train_confidence_function;
  
   
for endpointthreshold=.625:.125:1.125
      endtime=190;
    trackingparameters.endtime=endtime;
    endtime=endtime;
trackingparameters.endnodivthreshold=endpointthreshold;
tracking_driver_new_classifier_based_version;
if (~ROI)
    ROIxmin=0;
    ROIymin=0;
end
mkdir(outputdirectory);
saveGreedyNucleiFiles(esequence_trimmed,endtime,outputdirectory,anisotropy,ROIxmin,ROIymin);
    
uneddir=outputdirectory;

    test=evaluate_lineage_error(eddir,uneddir,endtime,anisotropy);
    test{37}=test{10}+test{21}+test{32};
    test{38}=test{12}+test{23}+test{34};
    test{39}=sum(test{37}+test{38});
    test{40}=sum(test{37});
    test{41}=sum(test{38});
    

errors{count}=test;
count=count+1
end