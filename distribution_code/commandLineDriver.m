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

%output
%outputs into embryodir/matlabnuclei a series of files with nuclear
%positions and diameters for each timepoint
%outputs into embryodir/embryonumber_fullmatlabresult.mat a dump of matlab
%detection results. The important part of this is the cell array esequence
%which contains structure fields finalpoints, finaldiams, which
%are the position and size of nuclei.

%closes matlab on completion since meant for batch mode remove exit below
%to prevent this


function commandLineDriver(parameterfile,embryodir,embryonumber,suffix)


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

xyres
anisotropy=zres/xyres*downsampling;
downsample=downsampling;
tlist=linspace(start_time,end_time,(end_time-start_time+1));
firsttimestep=start_time;
allvalid=[];


eall={};
tic
bottomdata={};

if(singlevolume&~nodata)
    zlevel=embryolevel;
else
    zlevel=slices;
end
%nucleibase=[nucleidir,embryonumber,'\'];
esequence={};
processSequence;
outputSNParamFile([embryodir,embryonumber,'_',suffix,'_sn_parameterfile.txt'],xyres,zres,slices,start_time,end_time,firsttimestepdiam)

toc


