%retraining procedure is to load a calculated mat result
%run this which regenerates the disk info data structure

% vis code run between comments here to bring nuclear image
% to use in correcting center list
%{
  circlevis=X*10;%Xorig;  
    circlevis(sub2ind(size(X),diskSet.xymax(:,2),diskSet.xymax(:,1),diskSet.xymax(:,3)))=-1;
    circlevis(sub2ind(size(X),filteredpointsswitched(:,2),filteredpointsswitched(:,1),filteredpointsswitched(:,3)))=1;


    for i=1:length(diskSet.centeredxymax)              
       % crustpoints=sphereFatCrustPoints_corrected(round(detdiameters(i)*1.2),filteredpointsswitched(i,:),anisotropy,Xorig,1);
        crustpoints=sphereBeltCrustPoints(xydetdiameters(i),diskSet.centeredxymax(i,:),anisotropy,X);
        circlevis(crustpoints)=-900;
    end
    
    for i=1:length(filteredpointsswitched)              
       % crustpoints=sphereFatCrustPoints_corrected(round(detdiameters(i)*1.2),filteredpointsswitched(i,:),anisotropy,Xorig,1);
        crustpoints=sphereBeltCrustPoints(detdiameters(i),filteredpointsswitched(i,:),anisotropy,X);
        circlevis(crustpoints)=2000;
    end
    
 for i=1:length(centers)


  points=centers{i};
     

    if (length(points)>0)
    points=diskSet.xymax(points(:,4),:);
    if length(points>0)
            points=points(:,1:3);
    circlevis(sub2ind(size(X),points(:,2),points(:,1),points(:,3)))=0;
    end
    end
 end
 nii=make_nii(circlevis)
 view_nii(nii)
 return
 
 
%}

%then does training
%{

%code for translating new data structures from pipeline to old variable
%names to run training code

%recalculate supp data for target timestep
time=277;%202
timepoint=time
previous=esequence{timepoint-1};
 if(newscope)
        if(~exist('LSM'))
         LSM=false;
        end

        if(LSM)
            X=imresize(loadCellStackLSM([embryodir,embryonumber],time,1,slices),downsample);

        else
            if rednuclei
                X=im2double((imresize(loadCellStackMetamorph([embryodir,embryonumber],time,1,slices,[0,0,0,0],zeropadding),downsample)));
                %end
                Xr=im2double((imresize(loadCellStackMetamorph([embryodir,embryonumber],time,2,slices,[0,0,0,0],zeropadding),downsample)));

            else
                X=im2double((imresize(loadCellStackMetamorph([embryodir,embryonumber],time,2,slices,[0,0,0,0],zeropadding),downsample)));
                %end
                Xr=im2double((imresize(loadCellStackMetamorph([embryodir,embryonumber],time,1,slices,[0,0,0,0],zeropadding),downsample)));
            end
        end
    else
        X=im2double(imresize(loadCellStack([embryodir,embryonumber,'_L1'],slices,time),downsample)); 
 end
Xorig=X;
 processVolume;
    X=tiledogfilter(Xorig,sigma,anisotropy);
 return
 %}    

centers=nucleiSet.centers;
%esequence{timepoint}.allcenters;

%bootstrap version without real hand labeled data
%{
range=nucleiSet.range;
corrected_centers=centers;
for i=1:length(centers)
    corrected_centers{i}=corrected_centers{i}(range{i},:);
end
good=linspace(1,length(centers),length(centers));
%}

%real version
corrected_centers=real_corrected_centers;
%good=[1,3,4,5,6,8,10,13,15,16,18,20,21,22,23,25,26,27,29,33,35,37,39,40,41,43,44,45,46,47,50,53,54,58,59,60,61,63,65,66,67,68,69,73,74,77,79,81,83,84,86,88,91,93,95,96,97,99,100,107,108,109,110,112,113,114,118,119,120,121,122,125,128,131,133,136,137,138,140,143,144,146,147,148,149,150,153,155,161,162,163,165,166,167,169,173,174,175,176,177,180,183,187,189,192,193,195,197,199,200,201,205,206,208,210,211,212,213,215,217,219,220,224,230,231,232,234,238,240,243,244,245,247,249,252,253,256,258,261,263,264];     
%hard ones started keeping track 35 ones actually had to change
%hardones=[37,45,53,54,58,61,68,73,77,81,84,86,88,91,93,95,96,97,99,100,107,112,113,118, 119,120,122,128,131,133,136,137,140,143,144,149,150,153,155,162,163,166,167,169,173,174,177,180,189,192,193,195,197,199,200,201,208,212,215,217,224,243,245,252,258,261,263,264];
%mr=[69,79]
%bad=[56,102,103];

%late stage good
good=[4,7,9,15,18,25,30,38,45,53,57,62,67,70,73,79,85,90,95,103,110,114,118,121,125,129,133,139,143,148,153,157,165,171,176,183,187,194,199,204,213,217,222,227,231,236,240,251,259,264,269,277,283,288,293,302,308,313,319,326,333,340,348,355,360,365,372,380,395,400,407,413,417,423,426,431,435,441,450,464,470];
good_ex=[4,6,7,9,10,15,16,18,19,25,26,30,31,38,39,45,46,53,54,57,58,62,63,67,68,70,71,73,74,79,85,86,90,91,95,96,103,110,111,114,,118,121,125,129,133,139,143,148,153,157,165,171,176,183,187,194,199,204,213,217,222,227,231,236,240,251,259,264,269,277,283,288,293,302,308,313,319,326,333,340,348,355,360,365,372,380,395,400,407,413,417,423,426,431,435,441,450,464,470];

oddones=[]; %ones that are not third on left but included anyway
bad=[380,475];
hardones=[30,53,57,73,90,103,125,129,133,139,153,157,165,171,199,204,213,227,231,236,240,251,264,269,277,283,288,293,302,308,313,319,340,348,355,365,372,390,395,400,413,417,423,435,457,464,470,477];
%80 is genuinely one disk 
filteredpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies,:);
originalpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies,:);

detdiameters=esequence{time}.diams;
xyzmaximavals=esequence{time}.maximavals;
xymaximavals=diskSet.xymaximavals;
xydetdiameters=diskSet.xydetdiameters;

% train distirbution of good and bad disks based on  a computed  result of
% disks centers and a hand corrected
%good2 is name indicating which center have been curated in
%corrected_centers
%old images good data
%good=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,34,35,36,37,38,39,40,41,42,43,44,45,46,47,162,163,164,165,166,167,168,169,170,171,172,174,175,176,177,178,179,181,182,183,184,185,186,187,237,238,239,240,241,242,244,245,247,248,249,250,251,252,253,254,255,256,257,258,259,260,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287;];

%slight doubt in my mind if this is general to whether center is included
%or not in set, but it should be


%figure
allgoodl=[];
allbadl=[]; %matricies of data for all features for all disks
allgoodr=[];
allbadr=[]; %matricies of data for all features for all disks
adjoining_bad_disks= 0;
hold on
oddcases=[];
minima_exists=0;
end_exists=0;
end_vals=[];
minima_vals=[];
not_minima=[];
for i=1:length(good)
    points=centers{good(i)};
    %** below comes out?
    points=[points(:,1:3);filteredpointsswitched(good(i),:)];

    zdiam=detdiameters(good(i));
    cpoints=corrected_centers{good(i)};
    points=centers{good(i)};
    sp=size(points);


    pivot=originalpointsswitched(good(i),3);
    topplane=min(min(points(:,3)),pivot);
    bottomplane=max(max(points(:,3)),pivot);

    %bounds on labeled data plus 3 in each direction is bound on
    %window to use as bad examples
    gtopplane=min(min(cpoints(:,3)),pivot)-2;
    gbottomplane=max(max(cpoints(:,3)),pivot)+2;
    %bounds on labeled data
    gtopplane=max(gtopplane,topplane);
    gbottomplane=min(gbottomplane,bottomplane);


    badones=[];
    badonesr=[];
    badonesl=[];
    goodonesr=[];
    goodonesl=[];

    counter=1;

    %calculate all values using live code
    features=calc_disk_feature_vector(points,xyzmaximavals(good(i)),detdiameters(good(i)),filteredpointsswitched(good(i),:),xymaximavals,xydetdiameters,anisotropy)   ;
    f2{good(i)}=features;


    for j=topplane:bottomplane-1
        leftdisks=points(find(points(:,3)==j),:);
        rightdisks=points(find(points(:,3)==j+1),:);
        if (j==pivot)
            %special case  left of pair is centerplane
 
            [m,minind]=min(distance(rightdisks(:,1:2)',filteredpointsswitched(good(i),1:2)'));
       
            %keep track of bad transitions to right rightdisk is bad
            if (length(find(rightdisks(minind,4)==cpoints(:,4)))==0 ) % a slice that was removed in correction
                %if(length(find(rightdisks(minind,4)==wronglyincluded(:,1)))>0) %fp
                if(rightdisks(minind,3)<=gbottomplane)
                    badones=[badones,counter];
                    badonesr=[badonesr,counter];
                end
            else
                goodonesr=[goodonesr,counter];
            end
            counter=counter+1;
        else
            if (j==pivot-1)
                % special case right of pair is centerplane
                sizes=size(leftdisks);

                for h=1:sizes(1)
                    if (length(find(leftdisks(h,4)==cpoints(:,4)))==0) % a slice that was removed in correction
                        %if(length(find(leftdisks(h,4)==wronglyincluded(:,1)))>0) %fp
                        if(leftdisks(h,3)>=gtopplane)
                            badones=[badones,counter];
                            badonesl=[badonesl,counter];
                        end
                    else
                        goodonesl=[goodonesl,counter];
                    end
                    counter=counter+1;
                end
            else
                %generic case
                sizes=size(leftdisks);
                for h=1:sizes(1)
                
                    if(j>pivot)
                        
                        [m,minind]=min(distance(rightdisks(:,1:2)',filteredpointsswitched(good(i),1:2)'));
       
                        if (length(find(rightdisks(minind,4)==cpoints(:,4)))~=0)
                            goodonesr=[goodonesr,counter];
                        end
                        if (length(find(rightdisks(minind,4)==cpoints(:,4)))==0)%&length(find(leftdisks(h,4)==cpoints(:,4)))~=0) % a slice that was removed in correction

                            if(rightdisks(minind,3)<=gbottomplane)
                                badones=[badones,counter];
                                badonesr=[badonesr,counter];
                            end
                        end
                    else
                        if (length(find(leftdisks(h,4)==cpoints(:,4)))~=0)
                            goodonesl=[goodonesl,counter];
                        end
                        if (length(find(leftdisks(h,4)==cpoints(:,4)))==0)%&length(find(rightdisks(minind,4)==cpoints(:,4)))~=0) % a slice that was removed in correction
                            if(leftdisks(h,3)>=gtopplane)
                                badones=[badones,counter];
                                badonesl=[badonesl,counter];
                            end
                        end
                    end
                    counter=counter+1;
                end
            end
        end

    end




    %assign featuers by different method
    allgoodl=[allgoodl,features(:,goodonesl)];
    allbadl=[allbadl,features(:,badonesl)];
    allgoodr=[allgoodr,features(:,goodonesr)];
    allbadr=[allbadr,features(:,badonesr)];



    points=centers{good(i)};
    points=[;filteredpointsswitched(good(i),:);points(:,1:3)];

    
end

  
  %scatter using 
  figure
  scatter3(allgoodl(4,:),allgoodl(5,:),allgoodl(6,:))
  hold on
    scatter3(allbadl(4,:),allbadl(5,:),allbadl(6,:),'r')
  
    scatter3(allgoodr(4,:),allgoodr(5,:),allgoodr(6,:),'b')

    scatter3(allbadr(4,:),allbadr(5,:),allbadr(6,:),'r')
    
        scatter3(allgoodlm(4),allgoodlm(5),allgoodlm(6),'c')
scatter3(allbadlm(4),allbadlm(5),allbadlm(6),'m')
 scatter3(allgoodrm(4),allgoodrm(5),allgoodrm(6),'c')
scatter3(allbadrm(4),allbadrm(5),allbadrm(6),'m')   
    
    
      figure
  scatter3(allgoodl(4,:),allgoodl(3,:),allgoodl(7,:))
  hold on
    scatter3(allbadl(4,:),allbadl(3,:),allbadl(7,:),'r')
  
    scatter3(allgoodr(4,:),allgoodr(3,:),allgoodr(7,:),'b')

    scatter3(allbadr(4,:),allbadr(3,:),allbadr(7,:),'r')
    
    scatter3(allgoodlm(4),allgoodlm(3),allgoodlm(7),'c')
scatter3(allbadlm(4),allbadlm(3),allbadlm(7),'m')
 scatter3(allgoodrm(4),allgoodrm(3),allgoodrm(7),'c')
scatter3(allbadrm(4),allbadrm(3),allbadrm(7),'m')   
    
          figure
  scatter3(allgoodl(4,:),allgoodl(1,:),allgoodl(2,:))
  hold on
    scatter3(allbadl(4,:),allbadl(1,:),allbadl(2,:),'r')
  
    scatter3(allgoodr(4,:),allgoodr(1,:),allgoodr(2,:),'b')

    scatter3(allbadr(4,:),allbadr(1,:),allbadr(2,:),'r')
    
scatter3(allgoodlm(4),allgoodlm(1),allgoodlm(2),'c')
scatter3(allbadlm(4),allbadlm(1),allbadlm(2),'m')
 scatter3(allgoodrm(4),allgoodrm(1),allgoodrm(2),'c')
scatter3(allbadrm(4),allbadrm(1),allbadrm(2),'m')   
    
%remove outliers of 9th round distribution
%{
allbadr=allbadr(:,find(allbadr(6,:)<1));
allbadl=allbadl(:,find(allbadl(6,:)>-1));
%}
    %remove outliers of 10th round distribution
    %{
allbadr=allbadr(:,find(allbadr(6,:)<.9&allbadr(5,:)>-.2));
allbadl=allbadl(:,find(allbadl(7,:)<.7));
allgoodr=allgoodr(:,find(allgoodr(1,:)<.4&allgoodr(5,:)<.2));
allgoodl=allgoodl(:,find(allgoodl(4,:)>-1));
%}
%now original training set

allbadr_u=allbadr;
allbadl_u=allbadl;
allgoodr_u=allgoodr;
allgoodl_u=allgoodl;



test=[1 2 3 4 5 6 7 ];
allbadr=allbadr_u(test,:);
allbadl=allbadl_u(test,:);
allgoodr=allgoodr_u(test,:);
allgoodl=allgoodl_u(test,:);

allbadrm=mean(allbadr');
allbadlm=mean(allbadl');
allgoodrm=mean(allgoodr');
allgoodlm=mean(allgoodl');
allgoodlc=cov(allgoodl');
allgoodrc=cov(allgoodr');
allbadlc=cov(allbadl');
allbadrc=cov(allbadr');

bbr=mvnpdf(allbadr',allbadrm,allbadrc);
bgr=mvnpdf(allbadr',allgoodrm,allgoodrc);

ggr=mvnpdf(allgoodr',allgoodrm,allgoodrc);
gbr=mvnpdf(allgoodr',allbadrm,allbadrc);

bbl=mvnpdf(allbadl',allbadlm,allbadlc);
bgl=mvnpdf(allbadl',allgoodlm,allgoodlc);

ggl=mvnpdf(allgoodl',allgoodlm,allgoodlc);
gbl=mvnpdf(allgoodl',allbadlm,allbadlc);

%correct left good
length(find(ggl>gbl))/length(ggl)
length(find(ggr>gbr))/length(ggr)
length(find(bbl>bgl))/length(bbl)
length(find(bbr>bgr))/length(bbr)
%wrong count
length(find(ggl<gbl))
length(find(ggr<gbr))
length(find(bbl<bgl))
length(find(bbr<bgr))

