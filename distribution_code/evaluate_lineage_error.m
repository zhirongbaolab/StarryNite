function errorInfo=evaluate_lineage_error(referencelineagedir,testlineagedir,endtime,anisotropy);
%generate error rates for detection and tracking given 2 directories
%containing 
numcells_labeled=[];

filtersulston=false;



%nucleidir='G:\My Documents\latetest\';


%embryonumbers = {'nuclei_matlab_besthack_2.5_noremove'};
%embryonumbers_c = {'edited_nuclei'};

allFN=[];
allFP=[];
downsample=1;
   numFP=[];
    numFN=[];
    nummatches=[];
    starttime=1;
   % endtime=280;%320;
    
 correct_sucessors=[];
 incorrect_sucessors=[];
numcello=[];
    
%xyres=.254;
%zres=1;   
%anisotropy=zres/xyres;

for time=starttime:endtime
 % for i=1:length(embryonumbers)

  %      embryonumber=embryonumbers{i};
  %      embryonumber_c=embryonumbers_c{i};
         nucleibase=testlineagedir;%[nucleidir,embryonumber,'\'];
       % nucleibase='G:\My Documents\canned_starrynite-1.3.3\src\nuclei\';
        nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei'];
        
        if~exist(nuclei,'file')
            'warning expected nuclei files missing'
            break
        end
       % [celldata,cellnames]=readnuclei_no_invalid(nuclei);
        [celldata_test1,cellnames]=readnuclei(nuclei);

        p1_sucessors_test=celldata_test1(:,9:10);
        celldiameter=mean(celldata_test1(:,7)); % size 
        matchthreshold=celldiameter*.5; 
        celllocations=celldata_test1(:,4:6);%pull nuclei from labeled data
        celllocations(:,1:2)=celllocations(:,1:2).*downsample; %compensate for downsampling
      
        
        
        % read corrected nuclei
       % nucleibase=[nucleidir,embryonumber,'\DOG_edited\nuclei\'];
         nucleibase=referencelineagedir;%[nucleidir,embryonumber_c,'\'];
        nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei']; 
        
        if~exist(nuclei,'file')
            'warning expected nuclei files missing'
            break
        end
        
        [celldata_c1,cellnames]=readnuclei(nuclei);
      
        %filter reference lineage to remove non sulston named nuclei  which
        %are sloppy editing
        if(filtersulston)
        goodednuclei=~(strncmp('Nuc',cellnames,3));
        celldata_c1=(celldata_c1(goodednuclei,:));
        end
        
        p1_sucessors=celldata_c1(:,9:10);
     
        s=size(celldata_c1);
        numcello=[numcello;s(1)];     
        celllocations_c=celldata_c1(:,4:6);%pull nuclei from labeled data
        
        %flip for comparing mirrored and unmirrored
        %celllocations_c(:,1)=celllocations_c(:,1)*-1+512;
        
        diameters=celldata_c1(:,7).*downsample;
        celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling
       
      
        
        [matches,matchessr]=compareDetectionWRadius_3_nonconflict(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);
  
        numcells_labeled=[numcells_labeled,length(diameters)];
        
        FP=celllocations(find(matches==-1),:);
        FN=celllocations_c(find(matchessr==-1),:);
        matchesc=celllocations(find(matches~=-1),:);   
        allFN=[allFN;time*ones(size(FN,1),1),FN];
        allFP=[allFP;time*ones(size(FP,1),1),FP];
        sizes=size(FP);
        numFP=[numFP,sizes(1)];
        sizes=size(FN);
        numFN=[numFN,sizes(1)];
        sizes=size(matchesc);
        nummatches=[nummatches,sizes(1)];
        
 %       cellspacing=...%mean(esequence{time}.selfdistances)-(mean(esequence{time}.finaldiameters)/2);
 %           mean(distance_anisotropic(celllocations_c',cellocations
    
%FNdetdensities=[FNdetdensities;cellspacing*ones(size(FN,1),1)];
%FPdetdensiteis=[FPdetdensities;cellspacing*ones(size(FP,1),1)];

        
    %tracking evaluation 
        if (time~=endtime)
           
            % read unedited nuclei
            %unedited SN

            %nucleibase=[nucleidir,embryonumber,'\DOG_unedited\nuclei\'];
                 %    nucleibase=[nucleidir,embryonumber,'\'];
                 nucleibase=testlineagedir;
                 nuclei=[nucleibase,'t',num2str(time+1,'%03d'),'-nuclei'];
                 
                 if~exist(nuclei,'file')
                     'warning expected nuclei files missing'
                     break
                 end
            [celldata_test2,cellnames]=readnuclei(nuclei);

            indiciesp2test=celldata_test2(:,2);
            celldiameter=mean(celldata_test2(:,7)); % size
            matchthreshold=celldiameter*.5;
            celllocations=celldata_test2(:,4:6);%pull nuclei from labeled data
            celllocations(:,1:2)=celllocations(:,1:2).*downsample; %compensate for downsampling


            % read corrected nuclei
            %nucleibase=[nucleidir,embryonumber,'\DOG_edited\nuclei\'];
             nucleibase=referencelineagedir;
               %nucleibase=[nucleidir,embryonumber_c,'\'];
            % nucleibase=[nucleidir,embryonumbers_c{i},'\'];
            nuclei=[nucleibase,'t',num2str(time+1,'%03d'),'-nuclei'];
            [celldata_c2,cellnames_c]=readnuclei(nuclei);
            
                    %filter reference lineage to remove non sulston named nuclei  which
        %are sloppy editing
        if(filtersulston)
        goodednuclei=~(strncmp('Nuc',cellnames,3));
        celldata_c1=(celldata_c1(goodednuclei,:));
        end
            
            indiciesp2=celldata_c2(:,2);
            celllocations_c=celldata_c2(:,4:6);%pull nuclei from labeled data
            diameters=celldata_c2(:,7).*downsample;
            celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling

            [matches2,matchessr2]=compareDetectionWRadius_3_nonconflict(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);

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
             [correct,wrong]=calculateHisFullTrackingError(p1_sucessors_t,p1_sucessors_test_t,matches,matchessr,matches2,matchessr2);
             
             correct_sucessors=[correct_sucessors;correct];
             incorrect_sucessors=[incorrect_sucessors;wrong];
             
           %  FNtrackdensities=[FNtrackdensities,celldensity*ones(sum(wrong(1:2)),1)];
           %  FPtrackdensities=[FPtrackdensities,celldensity*ones(sum(wrong(3:4)),1)];
        end
        
        % end
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
%   randind=round(rand(10,1)*size(allFP,1));
   
 %  allFP(round(rand(10,1)*size(allFP,1)),:)
   
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
      if(max(indearly)>length(correct_sucessors))
   indearly=indearly(indearly<=length(correct_sucessors));
   end
   
   results{10}=sum(incorrect_sucessors(indearly,1:4));
 %  results{11}=sum(incorrect_sucessors(indearly,1:4));
   %results{11}(1:2)= results{11}(1:2)./(sum(correct_sucessors(indearly,1)));
   %results{11}(3:4)= results{11}(3:4)./(sum(correct_sucessors(indearly,2)));
     results{11}(1:2)= (sum(correct_sucessors(indearly,1)));
   results{11}(3:4)= (sum(correct_sucessors(indearly,2)));
   
   results{12}=sum(incorrect_sucessors(indearly,5:8));
  % results{13}=sum(incorrect_sucessors(indearly,5:8));
  % results{13}(1:2)= results{13}(1:2)./(sum(correct_sucessors(indearly,3)));
  % results{13}(3:4)= results{13}(3:4)./(sum(correct_sucessors(indearly,4)));
     results{13}(1:2)=(sum(correct_sucessors(indearly,3)));
  results{13}(3:4)= (sum(correct_sucessors(indearly,4)));
     
     
   results{14}=(sum(numFN(indearly))+sum(nummatches(indearly)));
   
   results{15}='9th round error';
   results{16}=sum(numFN(ind9th))/(sum(numFN(ind9th))+sum(nummatches(ind9th)));
   results{17}=sum(numFP(ind9th))/(sum(numFN(ind9th))+sum(nummatches(ind9th)));
   results{18}=sum(numFN(ind9th)) ;
   results{19}=sum(numFP(ind9th));
   results{20}='tracing';
   
   %if this is end of list tracing has one less entry
   if(max(ind9th)>length(correct_sucessors))
   ind9th=ind9th(ind9th<=length(correct_sucessors));
   end
   
   results{21}=sum(incorrect_sucessors(ind9th,1:4),1);
%   results{22}=sum(incorrect_sucessors(ind9th,1:4),1);
%   results{22}(1:2)= results{22}(1:2)./(sum(correct_sucessors(ind9th,1),1));
%   results{22}(3:4)= results{22}(3:4)./(sum(correct_sucessors(ind9th,2),1));
   results{22}(1:2)= (sum(correct_sucessors(ind9th,1),1));
   results{22}(3:4)=(sum(correct_sucessors(ind9th,2),1));
   results{23}=sum(incorrect_sucessors(ind9th,5:8),1);
   %results{24}=sum(incorrect_sucessors(ind9th,5:8),1);
   %results{24}(1:2)= results{24}(1:2)./(sum(correct_sucessors(ind9th,3),1));
   %results{24}(3:4)= results{24}(3:4)./(sum(correct_sucessors(ind9th,4),1));
   
    results{24}(1:2)= (sum(correct_sucessors(ind9th,3),1));
   results{24}(3:4)= (sum(correct_sucessors(ind9th,4),1));
   
    
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
 results{32}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),1:4),1);
%  results{33}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),1:4),1);
 %  results{33}(1:2)= results{33}(1:2)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),1),1));
 %    results{33}(3:4)= results{33}(3:4)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),2),1));
    results{33}(1:2)=sum(correct_sucessors(ind10th(1:length(ind10th)-1),1),1);
     results{33}(3:4)= sum(correct_sucessors(ind10th(1:length(ind10th)-1),2),1);

    results{34}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),5:8),1);
 % results{35}=sum(incorrect_sucessors(ind10th(1:length(ind10th)-1),5:8),1);
 %  results{35}(1:2)= results{35}(1:2)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),3),1));
 %    results{35}(3:4)= results{35}(3:4)./(sum(correct_sucessors(ind10th(1:length(ind10th)-1),4),1));
   
   results{35}(1:2)= (sum(correct_sucessors(ind10th(1:length(ind10th)-1),3),1));
     results{35}(3:4)= (sum(correct_sucessors(ind10th(1:length(ind10th)-1),4),1));
   
   
   
results{36}=sum(numFN(ind10th))+sum(nummatches(ind10th));
results=results';


FNtrackdensities=[];
FPtrackdensities=[];
%generate count of errors of type at each timepoi

%nt can lookup density on
%other end
for i=1:size(incorrect_sucessors,1)

FNtrackdensities=[FNtrackdensities;(incorrect_sucessors(i,1:4))];
FPtrackdensities=[FPtrackdensities;(incorrect_sucessors(i,5:8))];
end
results{50}=numFN;
results{51}=numFP;
results{52}=FNtrackdensities;
results{53}=FPtrackdensities;
results{71}=numcells;
errorInfo=results;

   