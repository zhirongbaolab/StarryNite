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
    flipstack=true;
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
    
    
    if (MATLAB_STACK)
        load([embryodir,embryonumber,time_prefix,num2str(time,'%04d'),'.mat']);
        X=stack;%(:,:,1:slices);
        clear('stack');
        %X=X(:,:,1:slices);
        pause(1);%make easily ctrl-c able?
    else
        if(newscope)
            if rednuclei
                X=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,1,slices,[0,0,0,0],zeropadding))));
                %end
                Xr=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,2,slices,[0,0,0,0],zeropadding))));
                
            else
                X=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,2,slices,[0,0,0,0],zeropadding))));
                %end
                Xr=im2double(((loadCellStackMetamorph([embryodir,embryonumber],time,1,slices,[0,0,0,0],zeropadding))));
            end
            %our scope images are mirrored L,R
            
            if(flipstack)
                for s=1:size(X,3)
                    X(:,:,s)=fliplr(X(:,:,s));
                    Xr(:,:,s)=fliplr(Xr(:,:,s));
                end
            end
        else
            if(LSM_time)
                X=im2single(loadCellStackLSMtime([embryodir,embryonumber],time,LSM_channel,slices));
            else
                if (SIMPLETIFF)
                    name=[embryodir,embryonumber,time_prefix,num2str(time),'.tif'];
                    X=im2double(loadSimpleStackTiff(name));
                else
                    X=im2single((loadCellStack([embryodir,embryonumber],slices,time)));
                end
            end
        end
    end
    
    if(rednuclei)
        expind=2;
    else
        expind=1;
    end
    
    %output slices before subsampling via ROI
    if(outputSlice)
        %{
        if(slices>60||size(X,1)>600)%if big stack dont try to take percentile will run out of memory
            maximage=max(X,[],3);
            minval=prctile(reshape(maximage,[1,numel(maximage)]),outputSlice_linptilem);
            maxval=prctile(reshape(maximage,[1,numel(maximage)]),outputSlice_linptile);
            if(exist('Xr')&&tlist(example)==firsttimestep)
                
                Xrfinal=im2double(((loadCellStackMetamorph([embryodir,embryonumber],tlist(length(tlist)),expind,slices,[0,0,0,0],zeropadding))));
                maximage=max(Xrfinal,[],3);
                clear Xrfinal;
                minvalr=prctile(reshape(maximage,[1,numel(maximage)]),outputSlice_expptilem);
                maxvalr=prctile(reshape(maximage,[1,numel(maximage)]),outputSlice_expptile);
            end
        else
%}
    %        minval=prctile(reshape(X,[1,numel(X)]),outputSlice_linptilem);
    %        maxval=prctile(reshape(X,[1,numel(X)]),outputSlice_linptile);
            
            minval=approximateMatrixPercentile(X,outputSlice_linptilem,round((2^16)/10));
            maxval=approximateMatrixPercentile(X,outputSlice_linptile,round((2^16)/10));
   
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
  %      end
        %
        %if(rednuclei)
        %     outputAceTreeSlice(Xr,embryodir, embryonumber,time,minvalr,maxvalr,1,false);
        % else
        outputAceTreeSlice(X,embryodir, embryonumber,time,minval,maxval,1,false);
        %end
        if (exist('Xr'))
            %if(rednuclei)
            %    outputAceTreeSlice(X,embryodir, embryonumber,time,minval,maxval,1,true);
            %else
            outputAceTreeSlice(Xr,embryodir, embryonumber,time,minvalr,maxvalr,1,true);
            %end
            
        end
    end
    
    clear Xr;
    
    if(ROI)
        X=X(ROIymin:ROIymax,ROIxmin:ROIxmax,:);
    end
    
    X=imresize(X,downsample);
    
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

    backup_esequence=esequence;
    
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

if (SNoutput)
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
    
    if(ROI)
        outputSNFiles(embryodir,name,esequence,min(tlist),max(tlist),downsample,ROIxmin,ROIymin);
    else
         outputSNFiles(embryodir,name,esequence,min(tlist),max(tlist),downsample,1,1);

    end
        %else
    %    outputSNFiles([embryodir,'\',embryonumber],embryonumber,esequence,min(tlist),max(tlist),downsample);
    %end
end



if(exist('bigfile')&&bigfile==true)
    save([embryodir,name,'_fullmatlabresult.mat'],'-v7.3');
else
    save([embryodir,name,'_fullmatlabresult.mat']);
end
