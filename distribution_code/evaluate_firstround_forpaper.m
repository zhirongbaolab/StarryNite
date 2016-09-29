%generate starry nite error rates
numcells_labeled=[];




%nucleidir='H:\bao\dev\data analysis\data\';
%nucleidir='L:\bao\dev\data analysis\data\';


%nucleidir='G:\My Documents\posterdata\';


%embryonumbers = {'083105'};

    
%worm

  
 %mouse editing
nucleidir='L:\santella\mouse\Anthony\paperdata\test_data\';
embryonumbers = {'unedited'};
embryonumbers_c = {'edited'};
starttime=3;
endtime=5;
xyres=1.1;
zres=1.5;   

         %zebrafish 
nucleidir='L:\santella\zebrafish\paperdata\late_testdata\';
embryonumbers = {'unedited'};
embryonumbers_c = {'edited'};
starttime=401;%fake #s for acetree
endtime=403;   
xyres=1;
zres=5.18; 




  


%drosophila early
nucleidir='L:\santella\keller_d\paperdata\early_testdata\';
embryonumbers = {'unedited'};
embryonumbers_c = {'edited'};
xyres=.37;
zres=3.7; 
starttime=95;
    endtime=97;  
          starttime=44;
    endtime=46;
    
    
 nucleidir='G:\My Documents\latetest\';
embryonumbers = {'journalV_s1_besthacksettings'};
embryonumbers_c = {'journalV_s1_edit'};
starttime=20;
endtime=280
xyres=.254;
zres=1; 


downsample=1;
   numFP=[];
    numFN=[];
    nummatches=[];

    
 correct_sucessors=[];
 incorrect_sucessors=[];
numcello=[];
    
 allgaps_full={};
alldiams_full={};

anisotropy=zres/xyres;
alldiams=[];
alldiamss=[];
allgaps=[];
allgapss=[];
for time=starttime:endtime
  for i=1:length(embryonumbers)
      %{
           if (time>75)
         
        downsample=.5;
     else
            downsample=.25;
           end
     anisotropy=10*downsample;
      %}

      %{
        nucleibase=[nucleidir,embryonumbers{i},'\'];
        %read previous labeled timestep to simulate using previous step nuclear
        %diameter to set filter size
        nuclei=[nucleibase,'t',num2str(time-1,'%03d'),'-nuclei'];
        [celldata,cellnames]=readnuclei(nuclei);
        celldiameter=mean(celldata(:,7)); % size 
        matchthreshold=celldiameter*.5; 
        %}
       % slices=eachslice(i);
        embryonumber=embryonumbers{i};
                embryonumber_c=embryonumbers_c{i};
        % read unedited nuclei
       %unedited SN
      % nucleibase=[nucleidir,embryonumber,'\SN_unedited\nuclei\'];
      %      nucleibase=[nucleidir,embryonumber,'\DOG_unedited\nuclei\'];
        %%nucleibase=[nucleidir,embryonumber,'_unedited\'];
         nucleibase=[nucleidir,embryonumber,'\'];
       % nucleibase='G:\My Documents\canned_starrynite-1.3.3\src\nuclei\';
        nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei'];
       % [celldata,cellnames]=readnuclei_no_invalid(nuclei);
        [celldata_test1,cellnames]=readnuclei(nuclei);

        p1_sucessors_test=celldata_test1(:,9:10);
        celldiameter=mean(celldata_test1(:,7)); % size 
       
        matchthreshold=celldiameter*.5; 
        celllocations=celldata_test1(:,4:6);%pull nuclei from labeled data
        celllocations(:,1:2)=celllocations(:,1:2).*downsample; %compensate for downsampling
       gaps=mingaps(celllocations,celldata_test1(:,7),anisotropy);
              allgaps_full{time}=gaps;
              alldiams_full{time}=celldata_test1(:,7);
         alldiams=[alldiams;celldiameter];
          alldiamss=[alldiamss;std(celldata_test1(:,7))];
        allgaps=[allgaps;mean(gaps(gaps>0))];
       
      %{
      'using drosophila allgaps **'
      
      gaps=mingaps(celllocations,1.0634*celldata_test1(:,7),anisotropy);
      allgaps_full{time}=gaps;
        alldiams=[alldiams;mean(1.0634*celldata_test1(gaps'>=0&celllocations(:,3)>48&celllocations(:,1)>330,7))];
          alldiamss=[alldiamss;std(1.0634*celldata_test1(gaps'>=0&celllocations(:,3)>48&celllocations(:,1)>330,7))];
          
       allgaps=[allgaps;mean(gaps(gaps'>=0&celllocations(:,3)>48&celllocations(:,1)>330))];

       allgapss=[allgapss;std(gaps(gaps'>=0&celllocations(:,3)>48&celllocations(:,1)>330))];
 %}    
 %   allgaps=[allgaps;mean(gaps(gaps>=0))];
      % allgapss=[allgapss;std(gaps(gaps>=0))];

%{
        %note that lineaged version seems to be missing boundary nuclei
        
        %short circult and use first round data
       % celllocations=esequence{time-starttime+1}.firstroundpoints;
        celllocations=esequence{time-starttime+1}.finalpoints;
       
    
          bounds=[90,610,420,560,35,60];%fly expanded
               bounds=[90,610,400,580,35,60];%fly expanded2
                    
       %bounds=[0,460,170,260,25,45];% a mouse expanded
        bounds=[60,202,-1,202,0,12];%zebrafish late bounds
            bounds=[-1,402,-1,502,0,52];%zebrafish early bounds no cropping
             bounds=[90,610,400,590,35,60];%fly expanded3
       cshort=[];

        for j=1:length(celllocations)
    
            pcurrent=celllocations(j,:);
        %test bounds bc acetree editing seems broken on large data sets
            if(pcurrent(1)>bounds(1)&pcurrent(1)<bounds(2)&pcurrent(2)>bounds(3)&pcurrent(2)<bounds(4)&pcurrent(3)>bounds(5)&pcurrent(3)<bounds(6))
                cshort=[cshort;pcurrent];
            end
        end
 
        celllocations=cshort;
  %} 
         
        % read corrected nuclei
       % nucleibase=[nucleidir,embryonumber,'\DOG_edited\nuclei\'];
         nucleibase=[nucleidir,embryonumber_c,'\'];
        nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei']; 
        [celldata_c1,cellnames]=readnuclei(nuclei);
      
        p1_sucessors=celldata_c1(:,9:10);
     
        s=size(celldata_c1);
        numcello=[numcello;s(1)];     
        celllocations_c=celldata_c1(:,4:6);%pull nuclei from labeled data
        diameters=celldata_c1(:,7).*downsample;
        celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling
       
      
        
        [matches,matchessr]=compareDetectionWRadius_3(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);
  
        numcells_labeled=[numcells_labeled,length(diameters)];
        
        FP=celllocations(find(matches==-1),:);
        FN=celllocations_c(find(matchessr==-1),:);
        matches=celllocations(find(matches~=-1),:);   
        
        sizes=size(FP);
        numFP=[numFP,sizes(1)];
        sizes=size(FN);
        numFN=[numFN,sizes(1)];
        sizes=size(matches);
        nummatches=[nummatches,sizes(1)];
        
    
        if (time~=endtime)
           
            % read unedited nuclei
            %unedited SN

            %nucleibase=[nucleidir,embryonumber,'\DOG_unedited\nuclei\'];
                     nucleibase=[nucleidir,embryonumber,'\'];
             nuclei=[nucleibase,'t',num2str(time+1,'%03d'),'-nuclei'];
            [celldata_test2,cellnames]=readnuclei(nuclei);

            indiciesp2test=celldata_test2(:,2);
            celldiameter=mean(celldata_test2(:,7)); % size
            matchthreshold=celldiameter*.5;
            celllocations=celldata_test2(:,4:6);%pull nuclei from labeled data
            celllocations(:,1:2)=celllocations(:,1:2).*downsample; %compensate for downsampling


            % read corrected nuclei
            %nucleibase=[nucleidir,embryonumber,'\DOG_edited\nuclei\'];
               nucleibase=[nucleidir,embryonumber_c,'\'];
            % nucleibase=[nucleidir,embryonumbers_c{i},'\'];
            nuclei=[nucleibase,'t',num2str(time+1,'%03d'),'-nuclei'];
            [celldata_c2,cellnames]=readnuclei(nuclei);
            indiciesp2=celldata_c2(:,2);
            celllocations_c=celldata_c2(:,4:6);%pull nuclei from labeled data
            diameters=celldata_c2(:,7).*downsample;
            celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling

            [matches2,matchessr2]=compareDetectionWRadius_3(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);

             %translate sucessor indicies into array indicies (why is the
             %file formatted like this)
             s=size(p1_sucessors);
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
                 s=size(p1_sucessors_test);
             p1_sucessors_test_t=[];
             for j=1:s(1);
                 suc1=-1;
                 suc2=-1;
                 if(p1_sucessors_test(j,1)~=-1)
                     suc1=find(indiciesp2test==p1_sucessors_test(j,1));
                     if isempty(suc1)
                        suc1=-1;
                        'invalid nucleus pointed to by valid nucleus, should not happen in new files'
                     end
                 end
                  if(p1_sucessors_test(j,2)~=-1)
                     suc2=find(indiciesp2test==p1_sucessors_test(j,2));
                     if isempty(suc2)
                        suc2=-1;
                        'invalid nucleus pointed to by valid nucleus, should not happen in new files'
                     end
                 end
                 p1_sucessors_test_t=[p1_sucessors_test_t;suc1,suc2];
                 
                % p1_sucessors_test_t=[p1_sucessors_test_t;find(indiciesp2test==p1_sucessors_test(j,1)),find(indiciesp2test==p1_sucessors_test(j,2))];
             end
            
            [correct,wrong]=calculateTrackingError(p1_sucessors_t,p1_sucessors_test_t,matches,matchessr,matches2,matchessr2);

             correct_sucessors=[correct_sucessors;correct];
            incorrect_sucessors=[incorrect_sucessors;wrong];
        end
       
  end
end
 % hold on
 % plot(linspace(.4,2,3),ones(1,3)*mean(numFN./(numFN+nummatches)),'r')
 % plot(linspace(.4,2,3),ones(1,3)*mean(numFN./(numFN+nummatches))-std(numFN./(numFN+nummatches)),'r')  
 %  plot(linspace(.4,2,3),ones(1,3)*mean(numFN./(numFN+nummatches))+std(numFN./(numFN+nummatches)),'r') 
   
 %    plot(linspace(.4,2,3),ones(1,3)*mean(numFP./(numFN+nummatches)))
 % plot(linspace(.4,2,3),ones(1,3)*mean(numFP./(numFN+nummatches))-std(numFP./(numFN+nummatches)))  
  % plot(linspace(.4,2,3),ones(1,3)*mean(numFP./(numFN+nummatches))+std(numFP./(numFN+nummatches))) 
  %100 175
 
 %plot(175,mean(numFP./(numFN+nummatches)),'b*')
  % plot(175,mean(numFN./(numFN+nummatches)),'r*')
   
   %plot(mean(numFP./(numFN+nummatches)),mean(nummatches./(numFN+nummatches)),'g*')
   
   %note that numbers are mushed together by embryo, t*numebryos+embryoind
   %is location of particular result
   
   'global error'
   sum(numFN)/(sum(numFN)+sum(nummatches))   
   sum(numFP)/(sum(numFN)+sum(nummatches))
   %mean(numFP./(numFN+nummatches))
   
   numcells=numFN+nummatches;
   sum(numFN)
   sum(numFP)
   
   return
   indearly=find(numcells<194);
   ind9th=find(numcells>=194&numcells<=350);
   ind10th=find(numcells>350);
%ind10th=find(numcells>350&numcells<=422);
   
   'pre 9th round error'
   sum(numFN(indearly))/(sum(numFN(indearly))+sum(nummatches(indearly)))   
   sum(numFP(indearly))/(sum(numFN(indearly))+sum(nummatches(indearly)))
      sum(numFN(indearly))   
   sum(numFP(indearly))
    'tracing'
   sum(incorrect_sucessors(indearly))
   sum(incorrect_sucessors(indearly)/(sum(incorrect_sucessors(indearly))+sum(correct_sucessors(indearly))))
   
   (sum(numFN(indearly))+sum(nummatches(indearly)))
   
   '9th round error'
   sum(numFN(ind9th))/(sum(numFN(ind9th))+sum(nummatches(ind9th)))   
   sum(numFP(ind9th))/(sum(numFN(ind9th))+sum(nummatches(ind9th)))
   sum(numFN(ind9th)) 
   sum(numFP(ind9th))
   'tracing'
   sum(incorrect_sucessors(ind9th))
   sum(incorrect_sucessors(ind9th)/(sum(incorrect_sucessors(ind9th))+sum(correct_sucessors(ind9th))))
   
   (sum(numFN(ind9th))+sum(nummatches(ind9th)))
   '10th round error'
   sum(numFN(ind10th))/(sum(numFN(ind10th))+sum(nummatches(ind10th)))   
   sum(numFP(ind10th))/(sum(numFN(ind10th))+sum(nummatches(ind10th)))
   sum(numFN(ind10th))  
   sum(numFP(ind10th))
    'tracing'
   sum(incorrect_sucessors(ind10th(1:length(ind10th)-1)))
   sum(incorrect_sucessors(ind10th(1:length(ind10th)-1))/(sum(incorrect_sucessors(ind10th(1:length(ind10th)-1)))+sum(correct_sucessors(ind10th(1:length(ind10th)-1)))))
 
sum(numFN(ind10th))+sum(nummatches(ind10th))

   