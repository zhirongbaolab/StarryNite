%train confidence model
%match esequence contents agains corrected data
%and collect vector results
%note that though currently confidence is calculated from this model 
%it is never used the direct feature vector is, so this code mostly
%calculates the FP part of the answer key
filternames=false;%whether or not to filter out nuc named nuclei 


numcells_labeled=[];
%load('L:\duz\project\Imaging\BV24_NewJournalV_20091008\ZD_BV24_JournalV_1_s1_allinfo_logoddssum_fullmatlabresult.mat')
%nucleidir='G:\My Documents\latetest\';
%embryonumbers = {'nuclei_matlab_besthack_2.5_noremove'};
%embryonumbers_c = {'edited_nuclei'};
FPlocations={};

downsample=1;
FNFPdistances=[];
numFP=[];
numFN=[];
nummatches=[];
starttime=1;

%endtime=280;%135%280;%320;

%embryonumbers = {};
%nucleidir='L:\santella\unzipped_lineages\';
%embryonumbers_c={[lineages{lin},'_edited\nuclei\']};
%endtime=160;


correct_sucessors=[];
incorrect_sucessors=[];
numcello=[];

%xyres=.254;
%zres=1;
%anisotropy=zres/xyres;
global parameters;
parameters.anisotropyvector=[1,1,anisotropy];
FPdata=[];
TPdata=[];
for time=starttime:endtime
    
    %calculate confidence vectors
    confidencedata=calculateConfidenceVector(esequence{time},parameters);
    
    %embryonumber=embryonumbers{1};
    embryonumber_c=embryonumbers_c{1};
    
    
    celldiameter=mean(esequence{time}.finaldiams); % size
    matchthreshold=celldiameter*.5;
    celllocations=esequence{time}.finalpoints;%pull nuclei from labeled data
    
    
    if(ROI)
    celllocations(:,1)=celllocations(:,1)+ROIxmin;
    celllocations(:,2)=celllocations(:,2)+ROIymin;

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
    
    s=size(celldata_c1);
    celllocations_c=celldata_c1(:,4:6);%pull nuclei from labeled data
    diameters=celldata_c1(:,7).*downsample;
    celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling
    
    
    
    [matches,matchessr]=compareDetectionWRadius_3(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);
    esequence{time}.FP=(matches==-1);
    numcells_labeled=[numcells_labeled,length(diameters)];
    
    FP=celllocations(find(matches==-1),:);
    FN=celllocations_c(find(matchessr==-1),:);
    matchesc=celllocations(find(matches~=-1),:);
    
    FPdata=[FPdata;confidencedata(find(matches==-1),:)];
    TPdata=[TPdata;confidencedata(find(matches~=-1),:)];
    
    FPlocations{time}=esequence{time}.finalpoints(esequence{time}.FP,:);
    FNFPdistances=[FNFPdistances,min(distance(FP',FN'))./mean(diameters)];
end
%{

%%1 sumgfp, 2 nndistance (diams), 3 aspectratio, 4 logodds_avg, 5 avg gfp
figure
hold on
randone=round(rand([1,1000])*(size(TPdata,1)-1))+1;
scatter3(TPdata(randone,2),TPdata(randone,3),TPdata(randone,4),TPdata(randone,2)*10)
scatter3(FPdata(:,2),FPdata(:,4),FPdata(:,4),FPdata(:,2)*10)


figure
hold on
scatter(TPdata(randone,5),TPdata(randone,1))
scatter(FPdata(:,5),FPdata(:,1))

figure
hold on
randone=round(rand([1,1000])*(size(TPdata,1)-1))+1;
scatter3(TPdata(randone,2),TPdata(randone,3),TPdata(randone,1))
scatter3(FPdata(:,2),FPdata(:,3),FPdata(:,1))

%}
return

%confidence data values outside of this range are removed from training set
%to avoid outliers
%return
filtervalues=[inf,1.2,1.2,inf,8,2500];
goodvalues=min((FPdata<repmat(filtervalues,length(FPdata),1))');
FPdata=FPdata(goodvalues,:);
%filtervalues_good=[0,.5,0,0,0];
%goodgoodvalues=min((TPdata>repmat(filtervalues_good,length(TPdata),1))');
%TPdata=TPdata(goodgoodvalues,:);

con_goodmean=mean(TPdata);

con_goodstd=cov(TPdata);

con_badmean=mean(FPdata);
con_badstd=cov(FPdata);
trackingparameters.model.con_goodmean=con_goodmean;
trackingparameters.model.con_goodstd=con_goodstd;
trackingparameters.model.con_badmean=con_badmean;
trackingparameters.model.con_badstd=con_badstd;


%trained model
con_goodmean_1=mean(TPdata);

con_goodstd_1=cov(TPdata);



con_badmean_1=mean(FPdata);
con_badstd_1=cov(FPdata);

con_allmean_1=mean([FPdata;TPdata]);
con_allstd_1=cov([FPdata;TPdata]);
% end trained model
con_goodmean_2=mean(TPdata);

con_goodstd_2=cov(TPdata);



con_badmean_2=mean(FPdata);
con_badstd_2=cov(FPdata)




numFP=[];
numFN=[];
nummatches=[];
starttime=1%135;
endtime=280%280;%320;

correct_sucessors=[];
incorrect_sucessors=[];
numcello=[];

xyres=.254;
zres=1;
anisotropy=zres/xyres;
global parameters;
parameters.anisotropyvector=[1,1,anisotropy];
FPdata=[];
TPdata=[];
for time=starttime:endtime
    
    %calculate confidence vectors
    confidencedata=calculateConfidenceVector(esequence{time},parameters);
    
    embryonumber=embryonumbers{1};
    embryonumber_c=embryonumbers_c{1};
    
    
    celldiameter=mean(esequence{time}.finaldiams); % size
    matchthreshold=celldiameter*.5;
    celllocations=esequence{time}.finalpoints;%pull nuclei from labeled data
    
    
    
    % read corrected nuclei
    nucleibase=[nucleidir,embryonumber_c,'\'];
    nuclei=[nucleibase,'t',num2str(time,'%03d'),'-nuclei'];
    [celldata_c1,cellnames]=(readnuclei(nuclei));
    celldata_c1=double(celldata_c1);
    p1_sucessors=celldata_c1(:,9:10);
    
    s=size(celldata_c1);
    celllocations_c=celldata_c1(:,4:6);%pull nuclei from labeled data
    diameters=celldata_c1(:,7).*downsample;
    celllocations_c(:,1:2)=celllocations_c(:,1:2).*downsample; %compensate for downsampling
    
    
    
    [matches,matchessr]=compareDetectionWRadius_3(celllocations,celllocations_c,diameters*.5,1.5,anisotropy);
    esequence{time}.FP=(matches==-1);
    numcells_labeled=[numcells_labeled,length(diameters)];
    
    FP=celllocations(find(matches==-1),:);
    FN=celllocations_c(find(matchessr==-1),:);
    matchesc=celllocations(find(matches~=-1),:);
    
    FPdata=[FPdata;confidencedata(find(matches==-1),:)];
    TPdata=[TPdata;confidencedata(find(matches~=-1),:)];
end


filtervalues=[inf,1.2,1.2,inf,8,2500];
goodvalues=min((FPdata<repmat(filtervalues,length(FPdata),1))');
FPdata=FPdata(goodvalues,:);

con_goodmean_2=mean(TPdata);

con_goodstd_2=cov(TPdata);


con_badmean_2=mean(FPdata);
con_badstd_2=cov(FPdata);

con_allmean_2=mean([FPdata;TPdata]);
con_allstd_2=cov([FPdata;TPdata]);






trackingparameters.model.con_breakpoint=90;%number of cells at which switch btw models

trackingparameters.model.con_goodmean_1=con_goodmean_1;
trackingparameters.model.con_goodstd_1=con_goodstd_1;
trackingparameters.model.con_badmean_1=con_badmean_1;
trackingparameters.model.con_badstd_1=con_badstd_1;

trackingparameters.model.con_goodmean_2=con_goodmean_2;
trackingparameters.model.con_goodstd_2=con_goodstd_2;
trackingparameters.model.con_badmean_2=con_badmean_2;
trackingparameters.model.con_badstd_2=con_badstd_2;




%calculate confidences
for time=1:endtime
     confidencedata=calculateConfidenceVector(esequence{time},parameters);
  
     if(size(confidencedata,1)>trackingparameters.model.con_breakpoint)
        tempconfidences=mvnpdf(confidencedata, trackingparameters.model.con_goodmean_1,trackingparameters.model.con_goodstd_1)...
         ./mvnpdf(confidencedata, trackingparameters.model.con_badmean_1,trackingparameters.model.con_badstd_1);
     else
          tempconfidences=mvnpdf(confidencedata, trackingparameters.model.con_goodmean_2,trackingparameters.model.con_goodstd_1)...
         ./mvnpdf(confidencedata, trackingparameters.model.con_badmean_2,trackingparameters.model.con_badstd_2);
  
     end
 
   %  tempconfidences=confidencedata(:,1);
   %tempconfidences=mvnpdf(confidencedata, trackingparameters.model.con_badmean_1,trackingparameters.model.con_badstd_1);

   %mvnpdf(confidencedata, trackingparameters.model.con_goodmean_1,trackingparameters.model.con_goodstd_1);
   %...
    %     ./mvnpdf(confidencedata, trackingparameters.model.con_badmean_1,trackingparameters.model.con_badstd_1);

    tempconfidences=1./tempconfidences; %pdf ratio of bad
    %use raw confidences
    % tempconfidences(isnan(tempconfidences))=Inf;%...
     %   trackingparameters.model.confidence_confidencescaling_range;
     %tempconfidences(tempconfidences>trackingparameters.model.confidence_confidencescaling_range)=...
     %      trackingparameters.model.confidence_confidencescaling_range;
     %scale to bins
     %tempconfidences=tempconfidences./trackingparameters.model.confidence_confidencescaling_range*(trackingparameters.model.confidence_confidencescaling_buckets-1);
     %tempconfidences=trackingparameters.model.confidence_confidencescaling(round(tempconfidences)+1);
     if(max(isnan(tempconfidences)>0))
     'odd'
     end
     esequence{time}.confidences=tempconfidences;
     
end




return %

%over dims
figure
names={'log sum gfp','log nndist','aspect ratio','mean logodds claim','avg gfp'};
for i=1:5
    values=[TPdata(:,i);FPdata(:,i)];
his_buckets=[linspace(min(values),max(values),20)];
good_c=histc(TPdata(:,i),his_buckets);
bad_c=histc(FPdata(:,i),his_buckets);

subplot(1,5,i);
plot(his_buckets,good_c);
hold on
plot(his_buckets,bad_c,'r');
%plot(his_buckets,bad_c./(good_c+bad_c));

%axis([min(his_buckets),max(his_buckets),0,1])
title(names{i});

end


figure
for i=1:4
subplot(1,4,i)
hold on

scatter(TPdata(:,i),TPdata(:,5),'b');
scatter(FPdata(:,i),FPdata(:,5),'r');end



figure
for i=1:5
subplot(1,5,i)
hold on

scatter(TPdata(:,i),TPdata(:,6),'b');
scatter(FPdata(:,i),FPdata(:,6),'r');end


figure
hold on
scatter(TPdata(:,2),TPdata(:,3),'b');
scatter(FPdata(:,2),FPdata(:,3),'r');

num=linspace(1,length(FPdata),length(FPdata));
numtp=linspace(1,length(TPdata),length(TPdata));
figure
for i=1:5
subplot(1,5,i)
hold on

%scatter(TPdata(:,i),TPdata(:,6),num,'b');
scatter(FPdata(:,i),FPdata(:,6),num,'r');end


figure
hold on
scatter(TPdata(:,2),TPdata(:,3),numtp,'b');
scatter(FPdata(:,2),FPdata(:,3),num, 'r');

%{
bb=mvnpdf(FPdata,con_badmean,con_badstd);
bg=mvnpdf(FPdata,con_goodmean,con_goodstd);
ba=mvnpdf(FPdata,con_allmean,con_allstd);

gg=mvnpdf(TPdata,con_goodmean,con_goodstd);
gb=mvnpdf(TPdata,con_badmean,con_badstd);
ga=mvnpdf(TPdata,con_allmean,con_allstd);
%{
figure
scatter(gg,gb,'b');
hold on
scatter(bg,bb,'r');

figure
scatter(ga,gb,'b');
hold on
scatter(ba,bb,'r');


figure
scatter(ga,gg,'b');
hold on
scatter(ba,bg,'r');
%}
%correct left good
length(find(gg>gb))/length(gg)
length(find(bb>bg))/length(bb)
%wrong count
length(find(gg<gb))
length(find(bb<bg))

figure;hist(gg./ga);

figure;
scatter(gg,gb,'b');
hold on
scatter(bg,bb,'r');

ratio=length(gg)/length(bb); %how many times more good
%correct left good with prior weighting
length(find(gg*ratio>gb))/length(gg)
length(find(bb>bg*ratio))/length(bb)
%wrong count with prior weighting
length(find(gg*ratio<gb))
length(find(bb<bg*ratio))

goodhist=histc(gb./gg,linspace(0,7,25));
badhist=histc(bb./bg,linspace(0,7,25));
%figure
%hold on
%plot(linspace(0,10,100),goodhist,'b')
%plot(linspace(0,10,100),badhist,'r');



goodhist=histc(gb./gg,linspace(0,10,30));
badhist=histc(bb./bg,linspace(0,10,30));

goodpoints=(goodhist>2&badhist>2);
spacing=linspace(0,10,30)';
figure
%scatter(spacing(goodpoints),goodhist(goodpoints)./badhist(goodpoints),'r');
scatter(spacing(goodpoints),badhist(goodpoints)./goodhist(goodpoints),'r');

[estimates, model] = fitcurvedemo2(spacing(goodpoints),badhist(goodpoints)./goodhist(goodpoints))
[estimates, model] = fitcurvedemo2(spacing(goodpoints),goodhist(goodpoints)./badhist(goodpoints))



goodhist=histc(gb./gg,linspace(0,8,30));
badhist=histc(bb./bg,linspace(0,8,30));
goodhist(goodhist==0)=1;
badhist(badhist==0)=1;
result=badhist./goodhist;
result2=filter(ones(1,3)./3,1,result);





%vis of confidence
goodhist=histc(gb./gg,linspace(0,200,200));
badhist=histc(bb./bg,linspace(0,200,200));
d=[];
for i=length(goodhist):-1:1
    d=[d;sum(goodhist(i:end)),sum(badhist(i:end))];
end
figure
plot(d(:,1),d(:,2));
%}