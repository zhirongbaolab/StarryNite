%generate starry nite error rates
numcells_labeled=[];




%nucleidir='H:\bao\dev\data analysis\data\';
%nucleidir='L:\bao\dev\data analysis\data\';


%nucleidir='G:\My Documents\posterdata\';


%embryonumbers = {'083105'};



nucleidir='G:\My Documents\latetest\';


%embryonumbers = {'mirrored_wrong\journalV_s1_besthacksettings'};
%embryonumbers = {'mirrored_wrong\journalV_s1'};
%embryonumbers_c = {'mirrored_wrong\journalV_s1_edit'};

embryonumbers = {'nuclei_matlab_besthack_2.5_noremove'};
embryonumbers = {'tracktest_sn_nuclei'};
embryonumbers_c = {'edited_nuclei'};

%mouse editing
%nucleidir='L:\santella\mouse\Anthony\paperdata\training_data\';
%embryonumbers = {'unedited'};
%embryonumbers_c = {'edited'};


downsample=1;
   numFP=[];
    numFN=[];
    nummatches=[];
    starttime=1;
    endtime=280;%320;
    
 correct_sucessors=[];
 incorrect_sucessors=[];
numcello=[];
    
xyres=.254;
zres=1;   
anisotropy=zres/xyres;

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
      
        
        
        % read corrected nuclei
       % nucleibase=[nucleidir,embryonumber,'\DOG_edited\nuclei\'];
         nucleibase=[nucleidir,embryonumber_c,'\'];
        nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei']; 
        [celldata_c1,cellnames]=readnuclei(nuclei);
      
        p1_sucessors=celldata_c1(:,9:10);
     
        s=size(celldata_c1);
        numcello=[numcello;s(1)];     
        celllocations_c=celldata_c1(:,4:6);%pull nuclei from labeled data
        
        %flip for comparing mirrored and unmirrored
        %celllocations_c(:,1)=celllocations_c(:,1)*-1+512;
        
        diameters=celldata_c1(:,7).*downsample;
        celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling
       
      
        
        [matches,matchessr]=compareDetectionWRadius_3(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);
  
        numcells_labeled=[numcells_labeled,length(diameters)];
        
        FP=celllocations(find(matches==-1),:);
        FN=celllocations_c(find(matchessr==-1),:);
        matchesc=celllocations(find(matches~=-1),:);   
        
        sizes=size(FP);
        numFP=[numFP,sizes(1)];
        sizes=size(FN);
        numFN=[numFN,sizes(1)];
        sizes=size(matchesc);
        nummatches=[nummatches,sizes(1)];
        
    %tracking evaluation 
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
            
           % [correct,wrong]=calculateTrackingError(p1_sucessors_t,p1_sucessors_test_t,matches,matchessr,matches2,matchessr2);
[correct,wrong]=calculateFullTrackingError(p1_sucessors_t,p1_sucessors_test_t,matches,matchessr,matches2,matchessr2);

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
   results={};
   results{1}='global error';
   results{2}=sum(numFN)/(sum(numFN)+sum(nummatches));   
   results{3}=sum(numFP)/(sum(numFN)+sum(nummatches));
   %mean(numFP./(numFN+nummatches))
   
  numcells=numFN+nummatches;
   indearly=find(numcells<194);
   ind9th=find(numcells>=194&numcells<=350);
   ind10th=find(numcells>350);
%ind10th=find(numcells>350&numcells<=422);
   
   results{4}='pre 9th round error';
   results{5}=sum(numFN(indearly))/(sum(numFN(indearly))+sum(nummatches(indearly)))  ;
   results{6}=sum(numFP(indearly))/(sum(numFN(indearly))+sum(nummatches(indearly)));
   results{7}=   sum(numFN(indearly))   ;
   results{8}=sum(numFP(indearly));
   results{9}= 'tracing';
   results{10}=sum(incorrect_sucessors(indearly,1:4));
   results{11}=sum(incorrect_sucessors(indearly,1:4));
   results{11}(1:2)= results{11}(1:2)./(sum(correct_sucessors(indearly,1)));
   results{11}(3:4)= results{11}(3:4)./(sum(correct_sucessors(indearly,2)));
   results{12}=sum(incorrect_sucessors(indearly,5:8));
   results{13}=sum(incorrect_sucessors(indearly,5:8));
   results{13}(1:2)= results{13}(1:2)./(sum(correct_sucessors(indearly,3)));
   results{13}(3:4)= results{13}(3:4)./(sum(correct_sucessors(indearly,4)));
     
     
   results{14}=(sum(numFN(indearly))+sum(nummatches(indearly)));
   
   results{15}='9th round error';
   results{16}=sum(numFN(ind9th))/(sum(numFN(ind9th))+sum(nummatches(ind9th)))   ;
   results{17}=sum(numFP(ind9th))/(sum(numFN(ind9th))+sum(nummatches(ind9th)));
   results{18}=sum(numFN(ind9th)) ;
   results{19}=sum(numFP(ind9th));
   results{20}='tracing';
   results{21}=sum(incorrect_sucessors(ind9th,1:4));
   results{22}=sum(incorrect_sucessors(ind9th,1:4));
   results{22}(1:2)= results{22}(1:2)./(sum(correct_sucessors(ind9th,1)));
   results{22}(3:4)= results{22}(3:4)./(sum(correct_sucessors(ind9th,2)));
   results{23}=sum(incorrect_sucessors(ind9th,5:8));
   results{24}=sum(incorrect_sucessors(ind9th,5:8));
   results{24}(1:2)= results{24}(1:2)./(sum(correct_sucessors(ind9th,3)));
   results{24}(3:4)= results{24}(3:4)./(sum(correct_sucessors(ind9th,4)));
   
   
   %results{21}=sum(incorrect_sucessors(ind9th)/(sum(incorrect_sucessors(ind9th))+sum(correct_sucessors(ind9th))))
   
   results{25}=(sum(numFN(ind9th))+sum(nummatches(ind9th)));
   results{26}='10th round error';
   results{27}=sum(numFN(ind10th))/(sum(numFN(ind10th))+sum(nummatches(ind10th)))  ;
   results{28}=sum(numFP(ind10th))/(sum(numFN(ind10th))+sum(nummatches(ind10th)));
   results{29}=sum(numFN(ind10th))  ;
   results{30}=sum(numFP(ind10th));
   results{31}= 'tracing';
  % results{29}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1)))
  % results{32}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1))/(sum(incorrect_sucessors(ind10th(1:length(ind10th)-1)))+sum(correct_sucessors(ind10th(1:length(ind10th)-1)))))
 results{32}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),1:4));
  results{33}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),1:4));
   results{33}(1:2)= results{33}(1:2)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),1)));
     results{33}(3:4)= results{33}(3:4)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),2)));
    results{34}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),5:8));
  results{35}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),5:8));
   results{35}(1:2)= results{35}(1:2)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),3)));
     results{35}(3:4)= results{35}(3:4)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),4)));
   
   
   
results{36}=sum(numFN(ind10th))+sum(nummatches(ind10th));
results=results'

   