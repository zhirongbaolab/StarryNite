function [ esequence ] = matchLineageWithoutDeleted( esequence,trackingparameters,embryonumbers_c,nucleidir,ROI,ROIxmin,ROIymin )
%%iterate over computed lineage with corrected one in hand and compute
%%matching and correct successors 
filternames=false %flag to filter non sulston names in edited lineage, necessary in data sets which have unremoved cells from neighboring embyros or undeleted FP

downsample=1;
anisotropy=trackingparameters.anisotropyvector(3);
for time=trackingparameters.starttime:trackingparameters.endtime

    embryonumber_c=embryonumbers_c{1};
   
    celldiameter=mean(esequence{time}.finaldiams); % size
    matchthreshold=celldiameter*.5;
    celllocations=esequence{time}.finalpoints(~esequence{time}.delete,:);%pull nuclei from labeled data
    celllocations2=esequence{time+1}.finalpoints(~esequence{time+1}.delete,:);%pull nuclei from labeled data
    keeptpoints=~esequence{time}.delete;
    keept2points=~esequence{time+1}.delete;
    if(ROI)
    celllocations(:,1)=celllocations(:,1)+ROIxmin;
    celllocations(:,2)=celllocations(:,2)+ROIymin;
       celllocations2(:,1)=celllocations2(:,1)+ROIxmin;
    celllocations2(:,2)=celllocations2(:,2)+ROIymin;
    end
    
    % read corrected nuclei
    nucleibase=[nucleidir,embryonumber_c,'\'];
    nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei'];
    [celldata_c1,cellnames]=(readnuclei(nuclei));
    if (filternames)
                
                    %filter reference lineage to remove non sulston named nuclei  which
        %are sloppy editing
        goodednuclei=~(strncmp('Nuc',cellnames,3));
        celldata_c1=(celldata_c1(goodednuclei,:));
    end
    
    celldata_c1=double(celldata_c1);
    p1_sucessors=celldata_c1(:,9:10);
    
    
    indiciesp1=celldata_c1(:,2);
    s=size(celldata_c1);
    celllocations_c=celldata_c1(:,4:6);%pull nuclei from labeled data
    diameters=celldata_c1(:,7).*downsample;
    celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling
    
    [matches,matchessr]=compareDetectionWRadius_3_nonconflict(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);
    
    
    tmatches=ones(size(esequence{time}.finalpoints,1),1)*-1;
    tmatches(keeptpoints)=matches;
    esequence{time}.matches=tmatches;
    
    if (time<trackingparameters.endtime)
    %read t+1 nuclei
    nuclei=[nucleibase,'t',num2str(time+1,'%03d'),'-nuclei'];
    [celldata_c2,cellnames]=(readnuclei(nuclei));
    
       if (filternames)
                            
                    %filter reference lineage to remove non sulston named nuclei  which
        %are sloppy editing
        goodednuclei=~(strncmp('Nuc',cellnames,3));
        celldata_c2=(celldata_c2(goodednuclei,:));
       end    
    
    celldata_c2=double(celldata_c2);
    diameters2=celldata_c2(:,7).*downsample;
  
    celllocations_c2=celldata_c2(:,4:6);%pull nuclei from labeled data
    celllocations_c2(:,1:2)=celllocations_c2(:,1:2).*downsample; %compensate for downsampling
    indiciesp2=celldata_c2(:,2);
    [matches2,matchessr2]=compareDetectionWRadius_3_nonconflict(celllocations2,celllocations_c2,diameters2*.5,1.5,anisotropy);

    %translate successor ID indicies to row indicies
    p1_sucessors_t=[];
    for j=1:s(1);
        suc1=-1;
        suc2=-1;
        if(p1_sucessors(j,1)~=-1)
            suc1=find(indiciesp2==p1_sucessors(j,1));
            if isempty(suc1)
                suc1=-1;
                'invalid nucleus pointed to by valid nucleus, should not happen in new files'
            end
        end
        if(p1_sucessors(j,2)~=-1)
            suc2=find(indiciesp2==p1_sucessors(j,2));
            if isempty(suc2)
                suc2=-1;
                'invalid nucleus pointed to by valid nucleus, should not happen in new files'
            end
        end
        p1_sucessors_t=[p1_sucessors_t;suc1,suc2];
    end
    %store the translated successors for use in computing 'correct tracking answer'

    esequence{time}.corrected_sucessors=p1_sucessors_t;
    end
end

end

