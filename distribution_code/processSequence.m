if(~exist('NUC_FILE_OUTPUT'))
    NUC_FILE_OUTPUT=false;
end
if(~exist('LSM'))
    LSM=false;
end

if(~exist('LSM_channel'))
    LSM_channel=1;
end
if(~exist('SIMPLETIFF'))
    SIMPLETIFF=false;
end
if(~exist('LSM_time'))
    LSM_time=false;
end
if(~exist('MATLAB_STACK'))
    MATLAB_STACK=false;
end
if(~exist('flipstack'))
    flipstack=false;
end
if (~exist('outputSlice_linptile'))
    outputSlice_linptile=99.99;
    outputSlice_expptile=99;
    outputSlice_expptilem=25;
    outputSlice_linptilem=25;
end

for example=1:length(tlist)
    %put last volume data into place for use in local prediction
    if (tlist(example)==firsttimestep)
        previous=0;
    else
        previous=esequence{example-1};
    end

    time=tlist(example);
    e=[];
    if(~exist('time_prefix'))
        time_prefix='_';
    end
    
    %parse example file partially for matlab lsm and split cases
    %delegate this to matlab port of java name code for other cases so this is not going to be
    %used in simple tiff case except as an error correction last ditch if
    %the java name parsign doesn't generate an existing file name
    endprefixlocation=max(strfind(imageLocation,time_prefix));
    %if not found the parameter might be misconfigured try default
    if isempty(endprefixlocation)
        endprefixlocation=max(strfind(imageLocation,'_'));
    end
    if isempty(endprefixlocation)
        error('Can not identify time prefix string to parse image name aborting');
    end
    
    seplocation=max(strfind(imageLocation,filesep));
    embryodir=imageLocation(1:seplocation);
    embryonumber=imageLocation(seplocation+1:endprefixlocation-1);
    
    %note matlab stack is no longer explicit in new gui but mode will be tripped if 
    %selected file has .mat extension (old keller data)
    %and will keep this code back compatible with old gui should it
    %be needed
    if (MATLAB_STACK)   
        load([embryodir,embryonumber,time_prefix,num2str(time,'%04d'),'.mat']);
        X=stack;%(:,:,1:slices);
        clear('stack');
        %X=X(:,:,1:slices);
       % pause(1);%make easily ctrl-c able?
    else
        %old gui used newscope flag which was always a weird historical
        %choice, splitstack is used by new UI 
        %structured this way because newscope contaminates some old base
        %parameter files so it should only be paid attention to abscent a
        %more up to date directive. 
        if((~exist('splitstack','var')&&newscope)||splitstack)
            if rednuclei
                X=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,1,slices,[0,0,0,0],zeropadding))));
                %end
                Xr=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,2,slices,[0,0,0,0],zeropadding))));
                
            else
                X=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,2,slices,[0,0,0,0],zeropadding))));
                %end
                Xr=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,1,slices,[0,0,0,0],zeropadding))));
            end

        else
            %LSM no longer supported by new gui but left here for back
            %compatibility with new gui
            if(LSM_time)
                X=im2single(loadCellStackLSMtime([embryodir,embryonumber],time,LSM_channel,slices));
            else
                %now slice and implicit simpletiff are handled here
                %the current acetree.jar image handling code was ported
                %rather than called because it targets a newer java version
                %than is default in Matlab. 
                slicepattern='-p\d*.tif';
                isslice=~isempty(regexp(imageLocation,slicepattern));
                
                if(isslice)
                    %recompute prefix if slice bc sep is always -t 
                    endprefixlocation=max(strfind(imageLocation,'-t'));
                     embryonumber=imageLocation(seplocation+1:endprefixlocation-1);
   
                    %slice image
                    delim1='-';delim2='-';
                    X=loadCellStack([embryodir,embryonumber],slices,time,delim1,delim2);
                else
                    %simple tiff of isim dispim etc
                    %this call peels off the time
                    imagenameprefix=getImagePrefix(imageLocation);
    
                     if (zeropadding)
                        %name=[basename,'_t',num2str(time,'%03d'),'.TIF'];
                         imagefilename=[imagenameprefix,num2str(time,'%03d'),'.tif'];
                  
                    else
                        imagefilename=[imagenameprefix,num2str(time),'.tif'];
                  
                    end
                    
                    if(~exist(imagefilename,'file'))
                        %fallthrough backup if for some reason new acetree
                        %handling code wont try  the old suffix based tiff
                        %assumption as insurance against uncovered cases
                        imagefilename=[embryodir,embryonumber,time_prefix,num2str(time),'.tif'];
                    end
                      X=im2double(loadSimpleStackTiff(imagefilename));
                    if(slices<size(X,3))
                        X=X(:,:,1:slices);
                    end
                end
                %{
                if (SIMPLETIFF)
                    imagefilename=[embryodir,embryonumber,time_prefix,num2str(time),'.tif'];
                    X=im2double(loadSimpleStackTiff(imagefilename));
                    if(slices<size(X,3))
                        X=X(:,:,1:slices);
                    end
                else 
                    %slice case
                    delim1='_';delim2='-';
                    X=loadCellStack([embryodir,embryonumber],slices,time,delim1,delim2);
                    %X=im2single((loadCellStack([embryodir,embryonumber],slices,time)));
                end
                
                %}
            end
        end
    end
    
                %our scope images are mirrored L,R
            
            if(flipstack)
                for s=1:size(X,3)
                    X(:,:,s)=fliplr(X(:,:,s));
                    if (exist('Xr'))
                    Xr(:,:,s)=fliplr(Xr(:,:,s));
                    end
                end
            end
    
            
            
            
            %output slices before subsampling via ROI
            if(outputSlice)
                minval=approximateMatrixPercentile(imresize(X,downsample),outputSlice_linptilem,round((2^16)/10));
                maxval=approximateMatrixPercentile(imresize(X,downsample),outputSlice_linptile,round((2^16)/10));
                
        if(exist('splitstack','var')&&splitstack)
                    if(rednuclei)
                        expind=2;
                    else
                        expind=1;
                    end
                    
                    
                    %red mapping is calculated on the first volume and uniform
                    %afterward to allow quantitation of as opposed to nuclear channel
                    %mapped per volume to maximize visibility in the 8bit images
                    if(exist('Xr')&&tlist(example)==firsttimestep)
                        Xrfinal=im2double(((loadCellStackMetamorph([embryodir,embryonumber],tlist(length(tlist)),expind,slices,[0,0,0,0],zeropadding))));
                        
                        %minvalr=prctile(reshape(Xrfinal,[1,numel(Xrfinal)]),outputSlice_expptilem);
                        %maxvalr=prctile(reshape(Xrfinal,[1,numel(Xrfinal)]),outputSlice_expptile);
                        minvalr=approximateMatrixPercentile(Xrfinal,outputSlice_expptilem,round((2^16)/10));
                        maxvalr=approximateMatrixPercentile(Xrfinal,outputSlice_expptile,round((2^16)/10));
                        
                        
                        clear Xrfinal;
                    end
                end
                outputAceTreeSlice(X,embryodir, embryonumber,time,minval,maxval,1,false);
                
                if (exist('Xr'))
                    outputAceTreeSlice(Xr,embryodir, embryonumber,time,minvalr,maxvalr,1,true);
                end
            end
            
            %clear Xr;
            
    if(ROI)
        X=X(ROIymin:ROIymax,ROIxmin:ROIxmax,:);
        if(exist('Xr'))
        Xr=Xr(ROIymin:ROIymax,ROIxmin:ROIxmax,:);
        end
    end
    
  if downsample~=1
      X=imresize(X,downsample);
  end
    if(exist('Xr'))
        Xr=imresize(Xr,downsample);
    end
    
    svol=size(X);
    
    processVolume;
    
    
    %{
       % display code
        example=350;
        time=tlist(example);
        e=esequence{example};
        Xorig=im2double((imresize(loadCellStackMetamorph([embryodir,embryonumber],time,2,512,1024,30,[0,0,0,0]),downsample)));
       s=size(e.finalpoints);
        fakematch=linspace(1,s(1),s(1));
        displayData(fakematch,fakematch,round(e.finalpoints),round(e.finalpoints),Xorig)
    %}
    esequence{time}=e;
    
    ['processed ',num2str(example)]
    % clear unnecessary data files
    
    clear e;
    clear diskSet;
end

%bottom detection
if(~exist('bottomdetection'))
    bottomdetection=true;
end
if (bottomdetection&&~singlevolume&&nodatause)
    numcells=[];
    for i=1:length(esequence)
        if(~isempty(esequence{i}))
         numcells=[numcells;length(esequence{i}.finalpoints)];
        else
            numcells=[numcells;0];
        end
    end
    good=find((numcells>250)&(numcells<500));
    bottoms1=[];
    bottoms2=[];

       numonplane=zeros(1,svol(3));
    for i=1:length(good)
       % ind=find(esequence{good(i)}.finalmaximas>1.5*parameters.intensitythreshold);
       % planes=esequence{good(i)}.finalpoints(ind,3);
       % bottoms1=[bottoms1;max(planes)];
        planes=esequence{good(i)}.finalpoints(:,3);
     
         for j=1:svol(3) %ignore maxima on bottom 2 planes, they're artifacts
            numonplane(j)=numonplane(j)+length(find(planes==j));
         end
       
    end
    %bottom 2 definition is lowest plane with at least 20% of max number of
    %planes
    %modded this in attempt to replicate paper results...
    bottom=max(find(numonplane>=max(numonplane)*.2))+1;%with bottom maximas first small number plane seems to be real...
   
   % bottom=max(find(numonplane>=max(numonplane)*.2))+1;%with bottom maximas first small number plane seems to be real...
     %   bottom=max(find(numonplane>=max(numonplane)*.18));%can use this threshold if take them out
    %bottom=median(bottoms1);
    
%    bottomdata{emind}.bottom=bottom;
%    bottomdata{emind}.maxbottom=median(bottoms1);
%    bottomdata{emind}.bottoms=bottoms1;
%    bottomdata{emind}.numonplane=numonplane;

    clear X;
    clear Xr;

    % backup_esequence=esequence;
    
    for i=1:length(esequence)
        if(~isempty(esequence{i}))
            if(~isempty(esequence{i}.finalpoints))
                goodpoints=esequence{i}.finalpoints(:,3)<=bottom;
                esequence{i}.finalpoints=esequence{i}.finalpoints(goodpoints,:);
                esequence{i}.finaldiams=esequence{i}.finaldiams(goodpoints,:);
                esequence{i}.finalmaximas=esequence{i}.finalmaximas(goodpoints,:);
                esequence{i}.finalaveragepoints=esequence{i}.finalaveragepoints(goodpoints,:);
                esequence{i}.merged_sliceindicies=esequence{i}.merged_sliceindicies(goodpoints);
                esequence{i}.aspectratio=esequence{i}.aspectratio(goodpoints);
                
                    end
            
        end
        if (~nodata)
            %overwrite final with final after bottom detection
            [celldata,cellnames]=readnuclei([nucleidir,embryonumber_ed,'\t',num2str(tlist(i),'%03d'),'-nuclei']);
            threshold=.75*celldiameter; %for identity of detected and labeled nuclei
            celllocations=celldata(:,4:6);
            %Ucelllocations=Ucelldata(:,4:6);
            %stats from post bottom processing matching
            

            celllocations(:,1:2)=celllocations(:,1:2).*downsample;
            
            %match final merged and judged points
            filteredpointsswitched=esequence{i}.finalpoints;
            %compensate for ROI before matching against saved data set
            filteredpointsswitched(:,1)=filteredpointsswitched(:,1)+ROIxmin-1;
            filteredpointsswitched(:,2)=filteredpointsswitched(:,2)+ROIymin-1;
            [matches,matchessr]=compareDetectionWRadius_3(filteredpointsswitched,celllocations,celldata(:,7)*downsample*.5,1.5,anisotropy);
            
            FP=esequence{i}.finalpoints((matches==-1),:);
            FN=celllocations((matchessr==-1),:);
            esequence{i}.finalFN=FN;
            esequence{i}.finalFP=FP;
        end
    end
end

%if a polygonal ROI exists use it to filter the points in xy
%{
%move ROI 
if(exist('ROIpoints'))
    for i=1:length(esequence)
        if(~isempty(esequence{i}.finalpoints))
            pos=esequence{i}.finalpoints;
            polygon_points=ROIpoints;
            if(ROI)
                goodpoints=inpolygon(pos(:,1)+ROIxmin,pos(:,2)+ROIymin,polygon_points(:,1),polygon_points(:,2));
            else
                goodpoints=inpolygon(pos(:,1),pos(:,2),polygon_points(:,1),polygon_points(:,2));
           
            end
            esequence{i}.finalpoints=esequence{i}.finalpoints(goodpoints,:);
            esequence{i}.finaldiams=esequence{i}.finaldiams(goodpoints,:);
            esequence{i}.finalmaximas=esequence{i}.finalmaximas(goodpoints,:);
            esequence{i}.finalaveragepoints=esequence{i}.finalaveragepoints(goodpoints,:);
            esequence{i}.merged_sliceindicies=esequence{i}.merged_sliceindicies(goodpoints);
            esequence{i}.aspectratio=esequence{i}.aspectratio(goodpoints);

        end
    end
end
%}


%output results to SN 


%{
    alli=[];
    for i=1:length(esequence)
        alli=[alli;esequence{i}.finalmaximas];
    end

maxval=prctile(alli,99.99);
   maxval=prctile(alli,25); 
    
   
 
    
    
   
 
  esequence=backup_esequence;
    %}
    %hack for small diam bug in sn cant get smaller than 9 
    %{
smallest=[]
    
for i=1:length(esequence)
    esequence{i}.finaldiams=max(esequence{i}.finaldiams,9);
   smallest=[smallest,min(esequence{i}.finaldiams)];
end
%}
  %{
    
   %parameters.staging=[25,102,181,251,351];
   %parameters.intensitythreshold=[10,20,50,50,50,100];
    
    newthresh=100;
    xmax=421;
    ymax=315;
    zmax=26; %first illegal rather than last legal as above
    for i=1:length(esequence)
    numcells=length(esequence{i}.finalmaximas);
    if(numcells<25)
        newthresh=10;
    elseif(numcells<102) 
        newthresh=20;
    elseif(numcells<351)
        newthresh=75;
    else
        newthresh=125;
    end
    good=esequence{i}.finalmaximas>newthresh;
    good=good&esequence{i}.finalpoints(:,2)<xmax;
    good=good&esequence{i}.finalpoints(:,1)<ymax;
    good=good&esequence{i}.finalpoints(:,3)<zmax;
    esequence{i}.finaldiams=esequence{i}.finaldiams(good);
    esequence{i}.finalmaximas=esequence{i}.finalmaximas(good);
    esequence{i}.finalpoints=esequence{i}.finalpoints(good,:);
    end
     
    
  %}
    
    
    %if (newscope)
    if(isempty(suffix))
        name=embryonumber;
    else
       name=[embryonumber,'_',suffix];
    end
    %{
    %this is obsolete saving of positions for tracking by old C StarryNite
if (SNoutput)    
    if(ROI)
        outputSNFiles(embryodir,name,esequence,min(tlist),max(tlist),downsample,ROIxmin,ROIymin);
    else
         outputSNFiles(embryodir,name,esequence,min(tlist),max(tlist),downsample,1,1);

    end

end
    %}

