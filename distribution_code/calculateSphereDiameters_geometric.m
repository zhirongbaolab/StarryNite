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
            else
                if(~valley&&~negval&&val>minval) %valley
                    valley=offset-1;
                    allvalleydiameters(j,i)=offset-1;
                end
                
                
            end
            if(~negval&&(val<maxima*lowenough)) %gone near zero max-.8 max
                % if(~negval&&(val<maxima*.3)) %gone near zero max-.8 max
                negval=offset;
                alldiameters(j,i)=offset;%
                break
            end
        end
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
        [~,inval]=min(abs(alldistances(validprobes,i)-value));
        index=validprobes(inval); %index of the smallest valid probe
     
        for j=1:length(direction)-1 %for each ray except start which is assumed good
            nextindex=next(index);

            if (nextindex~=0)
                %if previous to judged is zero back previous index up till hit a nonzero value
                while(alldistances(index,i)==0)
                    index=previous(index);
                end
                rayratio=alldistances(nextindex,i)/alldistances(index,i);
                if (rayratio>largerthresh)||(rayratio<smallerthresh)
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

%edge bug that was never encountered bc of ROI
imagecoord_xpoints=min(max(1,imagecoord_xpoints),imsizes(2));
imagecoord_ypoints=min(max(1,imagecoord_ypoints),imsizes(1));


%[sortedrowsort,sortedrowsortindicies]=sort(raysort); %sorted version of ray position matrix
coverage=zeros(1,sizes(2));
centers_mod=zeros(sizes(2),3);
if(recenter)
    centers_mod=[];
    parfor i=1:sizes(2)
        coverage(i)=length(find(alldistances(:,i)~=0));
        %centroid
        %pull nonzero points and order them sucessively for centroid calculation
        %in practice this sorting doesnt acutally have significant time
        %contribution so I'm not going to fiddle with indexing to remove it
        orderindicies=raysort((alldistances(:,i)~=0)); %the relative location of the good indicies
        [~,sortindicies]=sort(orderindicies); %sorting these so rays go in right order
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
parfor i=1:sizes(2)
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


