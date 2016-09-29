%top level production driver for running lineage assuming that esequence
%and parameter data structure area already in current workspace

trackingparameters.trainingmode=false;
trackingparameters.recordanswers=false;


parameterConfiguration;
evalforced=false;

endtime=end_time;
trackingparameters.endtime=endtime;
trackingparameters.anisotropyvector=[1,1,anisotropy];
parameters.anisotropyvector=[1,1,anisotropy];

tracking_driver_new_classifier_based_version;
outputdirectory='./nuclei/';
mkdir(outputdirectory);
saveGreedyNucleiFiles(esequence,endtime,outputdirectory,anisotropy,ROIxmin,ROIymin);
zipname=[lineages{lin},'/',embryonumber,'_',suffix,'.zip'];
zip(zipname,[outputdirectory,'']);