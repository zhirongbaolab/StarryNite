
function [matches,correct_matched]=compareDetectionWRadius_3_nonconflict(detectedpoints,correctpoints,correctpointsizes,distancefactor,zscalingfactor)
%assigns matching labels to data based on radius
%returns detected long array of indicies into correctpoints
%a->b b->a less than *individual labeld radius

if isempty(detectedpoints)|isempty(correctpoints)

    matches=-1*ones(1,size(detectedpoints,1));
    correct_matched=-1*ones(1,size(correctpoints,1));
    return
end

%z values are in planes not pixels
s_detectedpoints=[detectedpoints(:,1),detectedpoints(:,2),detectedpoints(:,3)*zscalingfactor];
s_correctpoints=[correctpoints(:,1),correctpoints(:,2),correctpoints(:,3)*zscalingfactor];


distances=distance(s_detectedpoints',s_correctpoints');
if(size(detectedpoints,1)>1)
[mincorrect,Icorrect]=min(distances); %min distance and index of min dist for each correct point
else
  mincorrect=distances;
  Icorrect=ones(size(distances));
end
if(size(correctpoints,1)>1)
[mindetected,Idetected]=min(distances');%ditto for detected points
else
      mindetected=distances';
  Idetected=ones(size(distances'));
end
%I think this sorting is an artifact of previous version where matching was
%not bidirectional 

%sort mindetected,idetected so closest match gets first crack
[svals,sortIndex]=sortrows([mindetected',Idetected'],1);
smindetected=svals(:,1)'; %sorted list of min distance for each correct and index of detected
sIdetected=svals(:,2)';
     
 % 'modrule'

correct_matched=zeros(1,size(correctpoints,1));
%for each correct point 
FN=[];
matches=zeros(1,size(detectedpoints,1));
 %match a closest to b and distance ab is less than threhhold
 % because data is cleaned to 1 radius distance between cells
 %if threshold .5 radius this is stronger than match ab ba
 for i=1:length(smindetected)
     
     %detected point
     dpoint=s_detectedpoints(sortIndex(i),:);
     %correct point thats closes
     lpoint=s_correctpoints(sIdetected(i),:);
     
     %elliptical threshold distancefactor in z distancefactor/2 in x,y
   %  thresh_factor=correctpointsizes(sIdetected(i))*distancefactor;
    % thresh_val=(dpoint-lpoint).^2;
    % ovalratio=.75;
    % thresh_val=thresh_val(1)/(thresh_factor*ovalratio).^2+thresh_val(2)/(thresh_factor*ovalratio).^2+thresh_val(3)/(thresh_factor).^2;

    if( (smindetected(i)<correctpointsizes(sIdetected(i))*distancefactor) &&...               
        ( Icorrect(sIdetected(i))==sortIndex(i)|length(find(sIdetected==sIdetected(i)))==1))      
%Icorrect(sIdetected(i))==sortIndex(i));
        %is mutual match or at least no other corrected point wants this
        %detected point
        % 
     %now instead a-b b-a and less than elipse threshold
   %if(thresh_val<=1 && Icorrect(sIdetected(i))==sortIndex(i))  
        matches(sortIndex(i))=sIdetected(i);
        correct_matched(sIdetected(i))=sortIndex(i);
   else
       matches(sortIndex(i))=-1;
    end
end

for i=1:length(correct_matched)
    if(correct_matched(i)==0)
        correct_matched(i)=-1;
    end
    
end
 
 
 
%for i=1:length(mindetected)
   
%   if( (mindetected(i)<correctpointsizes(Idetected(i))*distancefactor) && correct_matched(Idetected(i))==0)  
%        matches(i)=Idetected(i);
%        correct_matched(Idetected(i))=i;
 %  else
 %      matches(i)=-1;
 %   end
%end

%for i=1:length(correct_matched)
%    if(correct_matched(i)==0)
%        correct_matched(i)=-1;
%    end
    
%end




