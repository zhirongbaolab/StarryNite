%runs nuclear detection 

%input
%parameterfile is full path location of parameter file
%imageLocation is an example image name, the one selected by user in ROI

%suffix is an optional string that will be appended to the name to all
%output keeping one processing of a set of image data separate from another
% for example to test multiple parameter sets, or (in conjunction with ROI parameters) 
%if multiple embryos are in one image 
%outputdirectory is where the output goes

%output
%outputs int outputdir/embryonumber_fullmatlabresult.mat a dump of matlab
%detection results. The important part of this is the cell array esequence
%which contains structure fields finalpoints, finaldiams, which
%are the position and size of nuclei.
%output int outputdir is a zip file contaning lineaged results

function detect_track_driver_allmatlab_v2(parameterfile,imageLocation,suffix,outputdirectory,polygon_points,isred);%,lineageParameterFileName)
runexpression=true;
newimage=false;

'beginning lineaging'
%cant profile compiled code
if (exist('profile'))
profile clear
profile on
end
tic
global parameters;


%cant pass blank string parameter on command line, so comes in as undef
if(~exist('suffix'))
    suffix='';
end


if(~exist('nodata'))
    nodata=true;%whether to use SN data to match online
end
nodatause=true;%whether to use SN data for diameter and stored bottom
savedata=true;
singlevolume=false;


readParameters;

load(distribution_file);

anisotropy=zres/xyres*downsampling;
downsample=downsampling;
tlist=linspace(start_time,end_time,(end_time-start_time+1));
firsttimestep=start_time;
allvalid=[];



eall={};

bottomdata={};

if(singlevolume&~nodata)
    zlevel=embryolevel;
else
    zlevel=slices;
end
%nucleibase=[nucleidir,embryonumber,'\'];
esequence={};
processSequence;

'detection completed, beginning lineaging'
set(0,'RecursionLimit',max(500,(end_time-start_time)*2));%if need more set recursion higher


mkdir([outputdirectory,suffix,embryonumber,'/nuclei']);


if (end_time-start_time>0)
    %parameterConfiguration;
    %parameterfile=lineageParameterFileName;
    %readParameters;
    trackingparameters.trainingmode=false;
    trackingparameters.recordanswers=false;
    evalforced=false;
    endtime=end_time;
    trackingparameters.endtime=endtime;
    trackingparameters.anisotropyvector=[1,1,anisotropy];
    parameters.anisotropyvector=[1,1,anisotropy];
    
    tracking_driver_new_classifier_based_version;

    %{
if(start_time~=1)
        tempesequence=esequence;
        esequence=cell(end_time,1);
        for i=1:length(tempesequence)
            esequence{start_time+i-1}=tempesequence{i};
        end
    end
%}    
    saveGreedyNucleiFiles(esequence,endtime,[outputdirectory,suffix,embryonumber,'/nuclei'],anisotropy,ROIxmin,ROIymin);
else
      %now done in detection 
      %{
    if(start_time~=1)
        tempesequence=esequence;
        esequence=cell(end_time,1);
        for i=1:length(tempesequence)
            esequence{start_time+i-1}=tempesequence{i};
        end
    end
      %}
    base=[outputdirectory,suffix,embryonumber,'/nuclei/'];
    output_unlineaged_acetree;
end

%move output here to save matlab version after tracking
%if(exist('bigfile')&&bigfile==true)
savematfile=true;
if savematfile
    save([outputdirectory,embryonumber,'_fullmatlabresult.mat'],'-v7.3');
end

%else
%    save([outputdirectory,embryonumber,'_fullmatlabresult.mat']);
%end

zipname=[outputdirectory,embryonumber,'_',suffix,'.zip'];
zipnameforfile=['./',embryonumber,'_',suffix,'.zip'];

zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%remove temp directory

%wrap because doesn't exist in compiled code.
if exist('profile')
profile viewer; 
end
rmdir([outputdirectory,suffix,embryonumber],'s');

toc

%then need to output xml file
%and run naming program
xmlname=[outputdirectory,embryonumber,'_',suffix,'.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');

%if (newimage)
    %fprintf (file,'<useStack type="1"/>\n');
   
    %fprintf (file,['<image file="',embryodir,embryonumber,'_t',num2str(start_time),'.tif"/>\n']);
%else
    %fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t',num2str(start_time,'%03d'),'-p01.tif"/>\n']);
%end
%note no longer need to explicitly distinguish slices etc it is auto detected
%but do need to specify if flipping or splitting in the interest of
%legiblity, though these will be typically interpreted correctly based on
%image name

cleanimagelocation=imageLocation;
cleanimagelocation(strfind(imageLocation,'\'))='/';
if(outputSlice)
    %if output slices save an xml file that opens them
    newDir=[embryodir,'image/tif'];
    newDir(strfind(newDir,'\'))='/';%avoid matlab escape charactering path
    cleanimagelocation=sprintf('%s/%s-t%03d-p%02d.tif',newDir,embryonumber,tlist(1),1);
    fprintf (file,['<image file="',cleanimagelocation,'"/>\n']);
    fprintf (file,'<Flip FlipMode="0"/>\n');
    fprintf (file,'<Split SplitMode="0"/>\n');
else
    fprintf (file,['<image file="',cleanimagelocation,'"/>\n']);
    
    if(~exist('flipstack')||~flipstack)
        fprintf (file,'<Flip FlipMode="0"/>\n');
    else
        fprintf (file,'<Flip FlipMode="1"/>\n');
    end
    if(~exist('splitstack')||~splitstack)
        fprintf (file,'<Split SplitMode="0"/>\n');
    else
        fprintf (file,'<Split SplitMode="1"/>\n');
    end
    
end

%if (time_prefix is not standard (ie has been manually set in param file
%propogate it to acetree
if(exist('time_prefix','var')&&~strcmp(time_prefix,'_'))
    fprintf (file,['<TimePrefix Prefix="',time_prefix,'"/>\n']);
end

fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,strcat('<end index="',num2str(end_time),'"/>\n'));
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);

%edited duplicate 
zipnameedited=[outputdirectory,embryonumber,'_',suffix,'_edited.zip'];
%ultimately this should be output but not doing it right now bc broken in 
%current test acetree
zipnameeditedforfile=['./',embryonumber,'_',suffix,'_edited.zip'];

xmlname=[outputdirectory,embryonumber,'_',suffix,'_edited.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');
%{
if (newimage)
    %if looking at new images delete old images after use
   % rmdir([embryodir,'image/'],'s');
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t1.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t001-p01.tif"/>\n']);
end
%}

%note no longer need to explicitly distinguish slices etc it is auto detected
%but do need to specify if flipping or splitting in the interest of
%legiblity, though these will be typically interpreted correctly based on
%image name
if(outputSlice)
    %if output slices save an xml file that opens them
    newDir=[embryodir,'image/tif'];
    newDir(strfind(newDir,'\'))='/';%avoid matlab escape charactering path
    cleanimagelocation=sprintf('%s/%s-t%03d-p%02d.tif',newDir,embryonumber,tlist(1),1);
    fprintf (file,['<image file="',cleanimagelocation,'"/>\n']);
    fprintf (file,'<Flip FlipMode="0"/>\n');
    fprintf (file,'<Split SplitMode="0"/>\n');
else
    
    fprintf (file,['<image file="',cleanimagelocation,'"/>\n']);
    if(~exist('flipstack')||~flipstack)
        fprintf (file,['<Flip FlipMode="0"/>\n']);
    else
        fprintf (file,['<Flip FlipMode="1"/>\n']);
    end
    if(~exist('splitstack')||~splitstack)
        fprintf (file,['<Split SplitMode="0"/>\n']);
    else
        fprintf (file,['<Split SplitMode="1"/>\n']);
    end
end

fprintf (file,['<nuclei file="',zipnameedited,'"/>\n']);

fprintf (file,strcat('<end index="',num2str(end_time),'"/>\n'));
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);
%edited duplicate 
copyfile (zipname,zipnameedited);

%currentdir=pwd;
%cd ('l:/bin/starryniteII');
%system([' java -Xmx500m -cp acebatch2.jar Measure1 ',xmlname]);
%2019 acebatch2 revisions
%call and surround xmlname with quotes in case contains spaces'
%actually took spaces out doesnt seem to work
system([' java -Xmx500m -jar acebatch2.jar Measure ',xmlname]);

%cd (currentdir);

'lineaging completed'

'running expression'
if(runexpression)
if(~isred)
     %system(['java -cp acebatch2.jar SixteenBitGreenExtractor1 ',xmlname,' 400']);  
  system(['java -jar acebatch2.jar Extractor ',xmlname,' G ',num2str(end_time)]);  
  
else
     %system(['java -cp acebatch2.jar SixteenBitRedExtractor1 ',xmlname,' 400']);  
     system(['java -cp acebatch2.jar Extractor ',xmlname,' R ',num2str(end_time)]);    
end
end
