function allbifurcationinfo= trainingDriver_functioncall(baseediteddir,matfilelineages,editedlineages,edittimes,newparameterfile)

%baseediteddir is directory just below that containign nuclie directory
%with edited nuclei
%matfilelineages is full path to mat file
%editedlineages the name of the directory containing edited lineages (to be
%combined with baseediteddir
%edittime last time edited in edited
%parameterfilenew optional parameter file containing settings to replace those in mat file 


%{
baseediteddir='L:\santella\blood_stemcell\usbkey\20120503\Lycorine_6_20120503_92523 AM\';
matfilelineages={'L:\santella\blood_stemcell\usbkey\20120503\Lycorine_6_20120503_92523 AM\Lycorine_6_emb_hysteresis_plusmodel1_fullmatlabresult.mat'}
editedlineages={'Lycorine_6_emb1_v3_edited'};
edittimes=189;
parameterfile='l:/santella/blood_stemcell/bloodstem_params.txt';
allbifurcationinfo= ...
trainingDriver_functioncall(baseediteddir,matfilelineages,editedlineages,edittimes,parameterfile)

%}

%baseediteddir is directory containing

%basematdir='L:\santella\blood_stemcell\usbkey\20120503\';
 %parameterfile='l:/santella/blood_stemcell/bloodstem_params.txt';
 %basedir=basematdir;
%edittimes=[189,189];


%lineages={'Lycorine_6_20120503_92523 AM\Lycorine_6_emb1_v3',...
%    'DMSO_11_20120503_93510 AM\DMSO_11_emb1'};
%matlineages={'Lycorine_6_emb_hysteresis_plusmodel1','DMSO_11_emb_hysteresis_plusmodel1'};

%lineagedir={'Lycorine_6_20120503_92523 AM','DMSO_11_20120503_93510 AM'}

%trainingmode refers to bif. classifier training not confidence
%uses oracle
trackingparameters.trainingmode=true;
trackingparameters.recordanswers=true;
answerkey=true;
outputtrimmed=false;
replacemodel=true;
allexpectedchange={};
allbifurcationinfo=[];

%this needs to be above because of ROI ie need to load it before
%load file and ROI that is used in file
%its needed only if parameter that effects training (like gap window size
%or easy tracking settings is different from that stored in initial
%sementation run
%if (~isempty(parameterfile))
%    readParameters;
%end
     
for lin=1:length(editedlineages)
    endtime=edittimes(lin);
  
    %'loader for files in zhuos directory'
     %load([basedir,lineagedir{lin},'/',lineages{lin},'_fullmatlabresult.mat']);
    'loader for training files'
     %load([basematdir,matlineages{lin},'_fullmatlabresult.mat']);

     load(matfilelineages{lin});

     %if want to ovverride any params in used result need to do it here
    %note that to not break things the parameterfile should not specify any
    %ROI info so as not to overrride the real ROI info just loaded with the
    %mat file. 
    if (~isempty(newparameterfile))
        parameterfile=newparameterfile;
    readParameters;
    'loaded paramfile'
    end
    evalforced=false;
     skipbifurcation=false;%needed to add this because tracking parameters in loaded file will now overwrite it
     trackingparameters.endtime=endtime;
    trackingparameters.anisotropyvector=[1,1,anisotropy];
    parameters.anisotropyvector=[1,1,anisotropy];
    
    embryonumbers = {};
    nucleidir=baseediteddir;
    embryonumbers_c={[editedlineages{lin},'\nuclei\']};

    if (answerkey)
        train_tracking_statistics_function;
        train_confidence_function;
    end
    endtime=endtime+5;
    trackingparameters.endtime=trackingparameters.endtime+5;
    
    trackingparameters.recordanswers=false;
    
    tracking_driver_new_classifier_based_version;
    
    %time adjustment back for final evaluation
    %endtime=endtime-9;
    %trackingparameters.endtime=trackingparameters.endtime-9;
    %reset back for final evaluation
    if(trackingparameters.trainingmode==false)
        endtime=endtime-6;
        trackingparameters.endtime=trackingparameters.endtime-6;
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
   
    if(trackingparameters.trainingmode)
        allexpectedchange{lin}=expected_corrections;
    end
      
end
   

           