%set diameter and possibly recenter based on geometry of zero crossing
function [diameters,centers,coverage,imagecoord_xpoints,imagecoord_ypoints]=calculateSphereDiameters_geometric(centers,Xfilt,expecteddiameter,recenter)

coverage=[];
%short circut to return if called with null parameters because bit trickly
%to rewrite actual processing correctly for that case
if(isempty(centers))
    diameters=[];
    coverage=[];
    imagecoord_xpoints=[];
    imagecoord_ypoints=[];
    return
end

direction=[0,1;
    0,-1;
    1,0;
    -1,0;
    0.7071,0.7071;
    -0.7071,-0.7071;
    0.7071,-0.7071;
    -0.7071,0.7071;
    0.8409,0.3483;
    -0.8409,-0.3483;
    -0.8409,0.3483;
    0.8409,-0.3483;
    0.3483,0.8409;
    -0.3483,-0.8409;
    0.3483, -0.8409;
    -0.3483,0.8409];

previous=[16,15,9,10,13,14,12,11,5,6,4,3,1,2,7,8]; %which direction index is previous to each if going clockwise
next=[13,14,12,11,9,10,15,16,3,4,8,7,5,6,2,1]; %which direction index is previous to each if going clockwise
raysort=[1,9,5,13,3,11,7,15,4,12,14,6,2,10,8,16];%location of each ray index in index sorted to go in clockwise direction

%altraysort=[1,16,7,12,3,9,5,13,2,15,8,11,4,10,6,14];

s=size(centers);
sizes=size(Xfilt);
%preallocate distances
alldiameters=zeros(length(direction),s(1));
allvalleydiameters=zeros(length(direction),s(1));

%lowenough=.25;
lowenough=getParameter('boundary_percent',length(centers));
%find valleys and zero crossings
for j=1:length(direction)
    for i=1:s(1)
        valley=0;
        center=centers(i,:);
        minval=Xfilt(center(2),center(1),center(3));
        negval=0;
        maxima=minval;
        %+offset in 2nd is right in x direction duh
        for offset=1:ceil(expecteddiameter/2*3)
            xpos=min(sizes(1),max(1,round(center(2)+offset*direction(j,1))));
            ypos=min(sizes(2),max(1,round(center(1)+offset*direction(j,2))));
            val=Xfilt(xpos,ypos,center(3));

            if(val<minval) %track min value
                minval=val;
            end

            if(~valley&&~negval&&val>minval) %valley
                valley=offset-1;
                allvalleydiameters(j,i)=offset-1;
            end

            if(~negval&&(val<maxima*lowenough)) %gone near zero max-.8 max
           % if(~negval&&(val<maxima*.3)) %gone near zero max-.8 max
                negval=offset;
                alldiameters(j,i)=offset;%
                break
            end

        end

        %if(~negval&&valley)
            %'set by valley';;
        %    alldiameters(j,i)=valley;%
        %end
    end %calculation of direction for each
end% end directions



%offsets in # of steps need to be turned into distances
sizes=size(alldiameters);
if (sizes(2)>0)
    alldistances=alldiameters.*repmat((direction(:,1).^2+direction(:,2).^2).^.5 ,1,sizes(2));
    allvalleydistances=allvalleydiameters.*repmat((direction(:,1).^2+direction(:,2).^2).^.5 ,1,sizes(2));
else
    alldistances=[];
end

%{

%check failure to find min/zero
length(find(alldistances==0))
[y,x]=ind2sub(size(alldistances),find(alldistances==0)');
[centers(x,:),y']


test=Xfilt;
test(sub2ind(size(Xfilt),centers(x,2),centers(x,1),centers(x,3)))=20;%1;
offsetcent=centers(x,:);
offsetcent(:,1)=offsetcent(:,1)+direction(y,2)*2;
offsetcent(:,2)=offsetcent(:,2)+direction(y,1)*2;
offsetcent=round(offsetcent);
offsetcent(:,1:2)=min(offsetcent(:,1:2),512);
offsetcent(:,1:2)=max(offsetcent(:,1:2),1);
test(sub2ind(size(Xfilt),offsetcent(:,2),offsetcent(:,1),offsetcent(:,3)))=-20;%.55;

test=Xfilt;
%draw all points
[y,x]=ind2sub(size(alldistances),find(alldistances~=0)');

%centers(:,1:2)=min(centers(:,1:2),512);
test(sub2ind(size(Xfilt),centers(x,2),centers(x,1),centers(x,3)))=100;%1;
offsetcent=centers(x,:);
for index=1:length(y) 
    offsetcent(index,1)=offsetcent(index,1)+direction(y(index),2)*alldistances(y(index),x(index))/(direction(y(index),1).^2+direction(y(index),2).^2)^.5;
    offsetcent(index,2)=offsetcent(index,2)+direction(y(index),1)*alldistances(y(index),x(index))/(direction(y(index),1).^2+direction(y(index),2).^2)^.5;
end
offsetcent=round(offsetcent);
offsetcent(:,1:2)=min(offsetcent(:,1:2),512);
offsetcent(:,1:2)=max(offsetcent(:,1:2),1);
test(sub2ind(size(Xfilt),offsetcent(:,2),offsetcent(:,1),offsetcent(:,3)))=-50;%-.55;


nii=make_nii(test)
view_nii(nii)
%}


%axisdiffs=[abs(alldistances(1,:)-alldistances(2,:));abs(alldistances(3,:)-alldistances(4,:));abs(alldistances(5,:)-alldistances(6,:));abs(alldistances(7,:)-alldistances(8,:))];
%axisratio=[abs(alldistances(1,:)./alldistances(2,:));abs(alldistances(3,:)./alldistances(4,:));abs(alldistances(5,:)./alldistances(6,:));abs(alldistances(7,:)./alldistances(8,:))];
%axisratio=min(axisratio,1./axisratio);

%fratio=[];
%for i=1:sizes(2)
%    fratio=[fratio,mean(axisratio(find(axisratio(:,i)~=0),i))];
%end

%detect outliers and mark them as invalid
%this is currently generic outlier detection more specific would be to
%compare each direction to 2 neibhoring direction vectors, which expect to
%be closer to the same even in asymmetric and off center conditions
%adjacentratios=[];

%largerthresh=3;%1.5;
%smallerthresh=1/3;
largerthresh=getParameter('large_ray_threshold',length(centers));
smallerthresh=getParameter('small_ray_threshold',length(centers));
outlierfilter=true;

if (outlierfilter)

    for i=1:sizes(2)%for each nucleus
        %min ray
        validprobes=find(alldistances(:,i)~=0); %nonzero indicies

        % [value,inval]=min(alldistances(validprobes,i)); %min of these
        value=median(alldistances(validprobes,i)); %min of these
        %closest value to returned as mean is index, as returns mean
        %of middle 2
        [value2,inval]=min(abs(alldistances(validprobes,i)-value));
        index=validprobes(inval); %index of the smallest valid probe
     
        for j=1:length(direction)-1 %for each ray except start which is assumed good
            nextindex=next(index);

            if (nextindex~=0)
                %if previous to judged is zero back previous index up till hit a nonzero value
                while(alldistances(index,i)==0)
                    index=previous(index);
                end
              
                if ((alldistances(nextindex,i)/alldistances(index,i)>largerthresh)||(alldistances(nextindex,i)/alldistances(index,i)<(smallerthresh)))
                    alldistances(nextindex,i)=0; %mark as invalid
                end
            end %if one to be judged is zero nothing to do
            index=nextindex;
        end
    end

end

%replace invalid with valley if present
alldistances(alldistances==0&allvalleydistances~=0)=allvalleydistances(alldistances==0&allvalleydistances~=0);

imsizes=size(Xfilt);

xpoints=alldistances.*repmat(direction(:,2),1,sizes(2));
ypoints=alldistances.*repmat(direction(:,1),1,sizes(2));
imagecoord_xpoints=alldistances.*repmat(direction(:,2),1,sizes(2))+repmat(centers(:,1)',length(direction),1);
imagecoord_ypoints=alldistances.*repmat(direction(:,1),1,sizes(2))+repmat(centers(:,2)',length(direction),1);

%[sortedrowsort,sortedrowsortindicies]=sort(raysort); %sorted version of ray position matrix
coverage=zeros(1,sizes(2));
centers_mod=zeros(sizes(2),3);
if(recenter)
    centers_mod=[];
    for i=1:sizes(2)
        coverage(i)=length(find(alldistances(:,i)~=0));
     
      
        %centroid
        %pull nonzero points and order them sucessively for centroid
        %calculation
        %in practice this sorting doesnt acutally have significant time
        %contribution so I'm not going to fiddle with indexing to remove it
        orderindicies=raysort((alldistances(:,i)~=0)); %the relative location of the good indicies
        [sorted,sortindicies]=sort(orderindicies); %sorting these so rays go in right order
        actualindicies=find(alldistances(:,i)~=0);
        actualindicies=actualindicies(sortindicies);
        
        xs=xpoints(actualindicies,i);
        ys=ypoints(actualindicies,i);
        area=0;
        cx=0;
        cy=0;
        %calculate centroid
        for j=1:length(xs)
            plus1=max(1,mod(j+1,length(ys)+1));
            area=area+xs(j)*ys(plus1)-ys(j)*xs(plus1);
            cx=cx+(xs(j)+xs(plus1))*(xs(j)*ys(plus1)-ys(j)*xs(plus1));
            cy=cy+(ys(j)+ys(plus1))*(xs(j)*ys(plus1)-ys(j)*xs(plus1));
        end
        area=area/2;
        cx=round(cx/(6*area));
        cy=round(cy/(6*area));
        %these are offsets because coordinate system is centered at old
        %center
        if (length(xs)<13)
            centers_mod(i,:)=centers(i,:);
        else
            centers_mod(i,:)=[centers(i,1)+cx,centers(i,2)+cy,centers(i,3)];
        end
    end
    centers_mod=min(max(1,centers_mod),repmat([imsizes(2),imsizes(1),imsizes(3)],sizes(2),1));

    %if not degenerate recenter otherwise dont risk changing it
    centers((coverage>13),:)=centers_mod((coverage>13),:);
end


%diameter=2*80th percentile of distances from new center
diameters=zeros(1,sizes(2));
for i=1:sizes(2)
    if(isempty(coverage)||coverage(i)>13)   
        
      %  time1=[];
      %  time2=[];
      %  for test=1:100
      %  tic
        %centerdistances=distance(centers(i,1:2)',[imagecoord_xpoints((alldistances(:,i)~=0),i),imagecoord_ypoints((alldistances(:,i)~=0),i)]');
     %   time1=[time1,toc];
     %    tic
        centerdistances=distanceToPoint(centers(i,1:2)',[imagecoord_xpoints((alldistances(:,i)~=0),i),imagecoord_ypoints((alldistances(:,i)~=0),i)]');
     %     time2=[time2,toc];
     %   end
        
        sdist=sort(centerdistances);
        diameters(i)=round(sdist(round(length(sdist)*.8)))*2+1;%80th percentile
    else
        diameters(i)=round(expecteddiameter);
    end
end


