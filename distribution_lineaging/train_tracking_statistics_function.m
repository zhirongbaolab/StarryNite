%train division nondivision model
%match esequence contents agains corrected data
%and collect vector results
filternames=false %whether or not to filter out nuc named nuclei 
lowdim=true %independent feature model instead of full cov
loosematch=false; %whether to use .5d or 2d a feature that isnt actualy necessary
%also cache matches and translated successor indicies in 
%esequence to be used in calculating correct computed nuclei centric 
%tracking answer
numcells_labeled=[];
%load('L:\duz\project\Imaging\BV24_NewJournalV_20091008\ZD_BV24_JournalV_1_s1_allinfo_logoddssum_fullmatlabresult.mat')
%nucleidir='G:\My Documents\latetest\';
%embryonumbers = {'nuclei_matlab_besthack_2.5_noremove'};
%embryonumbers_c = {'edited_nuclei'};
downsample=1;
starttime=1;
%endtime=280;%320;


%embryonumbers = {};
%nucleidir='L:\santella\unzipped_lineages\';
%embryonumbers_c={[lineages{lin},'_edited\nuclei\']};
%endtime=160;

%'nongeneral code running need to fix this !!!!!'
%xyres=.254;
%zres=1;
%anisotropy=zres/xyres;
interval=trackingparameters.interval;

global parameters;


%precompute per cell values
for time=starttime:endtime
    [integratedGFP,area]=integrateGFP(esequence{time},parameters);
    esequence{time}.totalGFP=integratedGFP;
    esequence{time}.avgGFP=integratedGFP./area;
end

divisions=0;
deaths=0;
realnuclei=0;
Divdata=[];
NoDivdata=[];
Tripledata=[];
for time=starttime:endtime
    
    
     distances=distance_anisotropic(esequence{time}.finalpoints',esequence{time}.finalpoints',trackingparameters.anisotropyvector);
 for j=1:max(size(distances))
    distances(j,j)=Inf;
 end
 mindistances=min(distances);
 esequence{time}.selfdistance=mindistances;
   
   % embryonumber=embryonumbers{1};
    embryonumber_c=embryonumbers_c{1};
   
    celldiameter=mean(esequence{time}.finaldiams); % size
    matchthreshold=celldiameter*.5;
    celllocations=esequence{time}.finalpoints;%pull nuclei from labeled data
    celllocations2=esequence{time+1}.finalpoints;%pull nuclei from labeled data

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
    
    if(filternames)
                
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
    
    dthresh=diameters*.5;
    if loosematch
        dthresh=diameters*2;
    end
    [matches,matchessr]=compareDetectionWRadius_3_nonconflict(celllocations,celllocations_c,...
        dthresh,1.5,anisotropy);
    
    esequence{time}.matches=matches;
    esequence{time}.matchesr=matchessr;
    realnuclei=realnuclei+length(p1_sucessors);
    if (time<endtime)
        deaths=deaths+length(find(p1_sucessors(:,1)==-1));
        divisions=divisions+length(find(p1_sucessors(:,2)~=-1));
    %read t+1 nuclei
    nuclei=[nucleibase,'t',num2str(time+1,'%03d'),'-nuclei'];
    [celldata_c2,cellnames]=(readnuclei(nuclei));
    
         if(filternames)       
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
    
     dthresh2=diameters2*.5;
    if loosematch
        dthresh2=diameters2*2;
    end
    
    [matches2,matchessr2]=compareDetectionWRadius_3_nonconflict(celllocations2,celllocations_c2,...
        dthresh2*.5,1.5,anisotropy);

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
    
    %compute and add to esequence gfp sum field
    [integratedGFP,area]=integrateGFP(esequence{time},parameters);
    esequence{time}.totalGFP=integratedGFP;
    esequence{time}.avgGFP=integratedGFP./area;
    
    %now walk through all tp and add entry to our training data set for
    %each  link which both ends of which have been successfuly detected
    anisotropy_vector=[1,1,anisotropy];
    for i=1:size(p1_sucessors,1) %iterate over real nuclei
     if(p1_sucessors_t(i,1)~=-1)
    if (matchessr(i)~=-1&&matchessr2(p1_sucessors_t(i,1))~=-1)%this cell and its sucessor was found
        %calculate pair vector for these two linked nuclei
        pairdata=calculateCellPairVector(esequence{time},matchessr(i),esequence{time+1},matchessr2(p1_sucessors_t(i,1)),anisotropy_vector);
	 pairdatanondivision=calculateCellPairVectorNondivision(esequence{time},matchessr(i),esequence{time+1},matchessr2(p1_sucessors_t(i,1)),anisotropy_vector);

        if(p1_sucessors_t(i,2)~=-1)
            Divdata=[Divdata;pairdata];
            if(matchessr2(p1_sucessors_t(i,2))~=-1) %second daughter also detected
                 pairdata=calculateCellPairVector(esequence{time},matchessr(i),esequence{time+1},matchessr2(p1_sucessors_t(i,2)),anisotropy_vector);
                  Divdata=[Divdata;pairdata];
                  Tripledata=[Tripledata;calculateCellTripleVector(esequence{time},matchessr(i),esequence{time+1},matchessr2(p1_sucessors_t(i,2)),esequence{time+1},matchessr2(p1_sucessors_t(i,1)),anisotropy_vector,interval)];
            end
        else
            NoDivdata=[NoDivdata;pairdatanondivision];           
        end
    
  
    end
     end
    end
    end %if not last timepoint
end

%end data collection
%vis of scatters



%trained model%%
div_mean=mean(Divdata);
div_std=cov(Divdata);
nodiv_mean=mean(NoDivdata);
nodiv_std=cov(NoDivdata);

div_triple_mean=mean(Tripledata);

if (lowdim)%learn independent feature instead of cov matrix
    div_triple_std=zeros(size(div_triple_mean,2));
    stddev=std(Tripledata);
    for ind=1:length(stddev)
        div_triple_std(ind,ind)=stddev(ind);
    end
else
    div_triple_std=cov(Tripledata);
end


esequence=computeTrackingAnswerKey(esequence,endtime);
if replacemodel

trackingparameters.model.div_mean=div_mean;
trackingparameters.model.div_std=div_std;
trackingparameters.model.nodiv_mean=nodiv_mean;
trackingparameters.model.nodiv_std=nodiv_std;
trackingparameters.model.div_triple_mean=div_triple_mean;
trackingparameters.model.div_triple_std=div_triple_std;

end
return


%for testing purposes copy correct answer to computed answer variables
for i=1:endtime
    esequence{i}.suc=esequence{i}.correct_suc;
    esequence{i}.suc_time=esequence{i}.correct_suc_time;
end
outputdirectory='l:\santella\nucleicorrectanswertest\nuclei\';
mkdir(outputdirectory);
saveGreedyNucleiFiles(esequence,endtime,outputdirectory,anisotropy);



figure
for i=2:4
subplot(1,3,i-1)
hold on
scatter(NoDivdata(:,i),NoDivdata(:,1),'b');
scatter(Divdata(:,i),Divdata(:,1),'g');
end

figure
for i=2:4
subplot(1,3,i-1)
hold on
scatter(Divdata(:,i),Divdata(:,1),'b');
end

figure
scatter3(NoDivdata(:,1),NoDivdata(:,2),NoDivdata(:,3),20*NoDivdata(:,4));
hold on
scatter3(Divdata(:,1),Divdata(:,2),Divdata(:,3),20*Divdata(:,4),'r');


figure
scatter3(NoDivdata(:,1),NoDivdata(:,2),NoDivdata(:,4));
hold on
scatter3(Divdata(:,1),Divdata(:,2),Divdata(:,4),'r');


%drift/elongation (1), velocity ratio(2), totalgfp ratio (3), avg %gfp(4),diameter(5)

tripledatanames={'drift/elongation','velocity ratio','total gfp ratio', 'avg gfp ratio', 'diameter'};
figure
v1=1;
v2=3;
v3=5;
scatter3(Tripledata(:,v1),Tripledata(:,v2),Tripledata(:,v3));
title('triple data');
xlabel(tripledatanames{v1});
ylabel(tripledatanames{v2});
zlabel(tripledatanames{v3});

figure
scatter3(Tripledata(:,1),Tripledata(:,2),Tripledata(:,4));
title('triple data');


%calculate the worst good score as a proxy value for when result is
%incalculable
divtriplescore=mvnpdf(Tripledata,div_triple_mean,div_triple_std);
divtriplescore=(1./divtriplescore);
values=log(divtriplescore);




divscore=mvnpdf(Divdata,div_mean,div_std);
divscore=double(1./divscore);
values=log(divscore);






nodivscore=mvnpdf(NoDivdata,nodiv_mean,nodiv_std);
nodivscore=(1./nodivscore);
values=log(nodivscore);

his_buckets=[-inf,linspace(prctile(values,1),prctile(values,85),50),inf];
his_buckets_display=[his_buckets(2)-1,linspace(prctile(values,1),prctile(values,85),50),his_buckets(51)+1];
good_c=histc(values,his_buckets);
figure
plot(his_buckets_display',good_c);

his_buckets=[-inf,linspace(prctile(values,1),prctile(values,95),20),inf];
his_buckets_display=[his_buckets(2)-1,linspace(prctile(values,1),prctile(values,95),20),his_buckets(21)+1];
good_c=histc(values,his_buckets);
figure
plot(his_buckets_display',good_c);




default_triple_score=1/prctile(divtriplescore,.05);


alt_div_score=1/prctile(divscore,.1);

%end trained model%
nodivscore=mvnpdf(NoDivdata,nodiv_mean,nodiv_std);

figure
hist(divscore)
title 'div pair score histogram'

figure
hist(nodivscore)
title 'non div score histogram';


