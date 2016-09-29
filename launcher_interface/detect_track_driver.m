%runs nuclear detection 

%input
%parameterfile is full path location of parameter file
% embryodir is full path of directory containing image files e.g.
% L:/disk1/bao/083105/image/tif/
% embryonumber is the prefix of the image file name e.g. 083105_L1
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

function detect_track_driver(parameterfile,embryodir,embryonumber,suffix,outputdirectory,polygon_points,isred)
newimage=false;%this flag controls whether we are outputting new images xml
% and deleting old images

'beginning lineaging'

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

%output sn param file
paramfile=[outputdirectory,embryonumber,'_',suffix,'_sn_parameterfile.txt'];
outputSNParamFile(paramfile,xyres,zres,slices,start_time,end_time,firsttimestepdiam)
nuclocation=[embryodir,embryonumber,'_',suffix,'_matlabnuclei\'];

mkdir([outputdirectory,suffix,embryonumber,'/nuclei']);

%lineage processing
%run starrynite

imagefile=[embryodir,'image/tif/',embryonumber];
['running ',snpath,' ',imagefile,' ',[outputdirectory,suffix,embryonumber,'/nuclei/'],'  ',paramfile,' ',nuclocation]
currentdir=pwd;
cd ([outputdirectory,suffix,embryonumber,'/nuclei/']);
system([snpath,' ',imagefile,' ',[outputdirectory,suffix,embryonumber,'/nuclei/'],'  ',paramfile,' ',nuclocation]);
%cd ([outputdirectory,suffix,embryonumber]);
cd '../../';
zipname=[outputdirectory,embryonumber,'_',suffix,'.zip'];
zip(zipname,[outputdirectory,suffix,embryonumber,'/nuclei']);
%remove temp directory

rmdir([outputdirectory,suffix,embryonumber],'s');







%then need to output xml file
%and run naming program
xmlname=[outputdirectory,embryonumber,'_',suffix,'.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');

if (newimage)
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t1.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t001-p01.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);

%edited duplicate 
zipnameedited=[outputdirectory,embryonumber,'_',suffix,'_edited.zip'];
xmlname=[outputdirectory,embryonumber,'_',suffix,'_edited.xml'];
file=fopen(xmlname,'w');

fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');
if (newimage)
    %if looking at new images delete old images after use
  %  rmdir([embryodir,'image/'],'s');
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t1.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t001-p01.tif"/>\n']);
end

fprintf (file,['<nuclei file="',zipnameedited,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);
%edited duplicate 
copyfile (zipname,zipnameedited);
cd(currentdir);%return to where were hoping acebatch2 is in path
%currentdir=pwd;
%cd ('l:/bin/starryniteII');
system([' java -Xmx500m -cp acebatch2.jar Measure1 ',xmlname]);
cd (currentdir);

'lineaging completed'
'running expression'
if(isred)
   system(['java -cp acebatch2.jar SixteenBitGreenExtractor1 ',xmlname,' 400']);
else
     
        system(['java -cp acebatch2.jar SixteenBitRedExtractor1 ',xmlname,' 400']);  
end

