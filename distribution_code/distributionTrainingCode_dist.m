%retraining procedure is to load a calculated mat result
%run this which regenerates the disk info data structure


%{
time=277;

%RUN THIS FIRST
%code for translating new data structures from pipeline to old variable
%names to run training code

%recalculate supp data for target timestep

timepoint=time;
previous=esequence{timepoint-1};

%load the image for the timepoint to be processed
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

filteredpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies,:);
originalpointsswitched=diskSet.centeredxymax(nucleiSet.centerindicies,:);

detdiameters=diskSet.xydetdiameters(nucleiSet.centerindicies);%esequence{time}.diams;
xyzmaximavals=diskSet.xymaximavals(nucleiSet.centerindicies);%esequence{time}.maximavals;
xymaximavals=diskSet.xymaximavals;
xydetdiameters=diskSet.xydetdiameters;


% vis code to bring up nuclear images
% to use in correcting center list
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
 centers=nucleiSet.centers;
corrected_centers=centers;    
    
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



% good corrected nuclei
good=[];


% train distirbution of good and bad disks based on  a computed  result of
% disks centers and a hand corrected

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

