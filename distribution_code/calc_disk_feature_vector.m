function features=calc_disk_feature_vector(points,intensity,diameter,centerpoint,xymaximavals,xydetdiameters,anisotropy)   
%outputs feature vector in disk order 

pivot=centerpoint(3);
topplane=min(min(points(:,3)),pivot);
bottomplane=max(max(points(:,3)),pivot);
s=size(points);

goodmaximasp=zeros(1,s(1));
gooddiamsp=zeros(1,s(1));
c_xydifp=zeros(1,s(1));
goodmaximadifsp=zeros(1,s(1));
goodmaximadifsp2=zeros(1,s(1));
xydisk_distancediff=zeros(1,s(1));
newzdist=zeros(1,s(1));

counter=1;
for j=topplane:bottomplane-1
    leftdisks=points((points(:,3)==j),:);
    rightdisks=points((points(:,3)==j+1),:);
    sizes=size(leftdisks);
    sizes2=size(rightdisks);

    if (j==pivot)
        %special case  left of pair is centerplane
        %leftdisks set containing the set of center disks if not just
        %centerpoint is ignored
        for h=1:sizes2(1)
            goodmaximasp(counter)=xymaximavals(rightdisks(h,4))./intensity;
            gooddiamsp(counter)=xydetdiameters(rightdisks(h,4))/diameter;
           % c_xydifp(counter)=distance(rightdisks(h,1:2)',centerpoint(1:2)')./diameter;
            
            c_xydifp(counter)=(rightdisks(h,1)-centerpoint(1)).^2;
            c_xydifp(counter)=c_xydifp(counter)+(rightdisks(h,2)-centerpoint(2)).^2;
            c_xydifp(counter)=(c_xydifp(counter)^.5)/diameter;
            
            goodmaximadifsp(counter)=(xydetdiameters(rightdisks(h,4))-diameter)./diameter;
            goodmaximadifsp2(counter)=(xymaximavals(rightdisks(h,4))-intensity)/intensity;
            
            xydisk_distancediff(counter)=(rightdisks(h,1)-centerpoint(1)).^2;
            xydisk_distancediff(counter)=xydisk_distancediff(counter)+(rightdisks(h,2)-centerpoint(2)).^2;
            xydisk_distancediff(counter)=(xydisk_distancediff(counter)^.5)/diameter;
            
           % xydisk_distancediff(counter)=(distance(rightdisks(h,1:2)',centerpoint(1:2)')-0)./(diameter);
            newzdist(counter)=(j+1-pivot)*anisotropy/(diameter);
            counter=counter+1;
        end
    else
        if (j==pivot-1)
            % special case right of pair is centerplane
            %rightdisks set containing the set of center disks if not just
            %centerpoint is ignored
            for h=1:sizes(1)
                goodmaximasp(counter)=xymaximavals(leftdisks(h,4))./intensity;
                gooddiamsp(counter)=xydetdiameters(leftdisks(h,4))/diameter;
                %c_xydifp(counter)=distance(leftdisks(h,1:2)',centerpoint(1:2)')./diameter;
                
                c_xydifp(counter)=(leftdisks(h,1)-centerpoint(1)).^2;
                c_xydifp(counter)=c_xydifp(counter)+(leftdisks(h,2)-centerpoint(2)).^2;
                c_xydifp(counter)=(c_xydifp(counter)^.5)/diameter;


                goodmaximadifsp(counter)=(diameter-xydetdiameters(leftdisks(h,4)))./diameter;
                goodmaximadifsp2(counter)=(intensity-xymaximavals(leftdisks(h,4)))/intensity;
                %one term cuts out bc distance is always zero
               % xydisk_distancediff(counter)=distance(leftdisks(h,1:2)',centerpoint(1:2)')./(diameter);
               
               xydisk_distancediff(counter)=(centerpoint(1)-leftdisks(h,1)).^2;
               xydisk_distancediff(counter)=xydisk_distancediff(counter)+(centerpoint(2)-leftdisks(h,2)).^2;
               xydisk_distancediff(counter)=(xydisk_distancediff(counter)^.5)/diameter;


                newzdist(counter)=(j-pivot)*anisotropy/(diameter);
                counter=counter+1;
            end
        else
            %generic case
            sizes=size(leftdisks);
            if(j<pivot)
                for h=1:sizes(1)
                    goodmaximasp(counter)=xymaximavals(leftdisks(h,4))./intensity;
                    gooddiamsp(counter)=xydetdiameters(leftdisks(h,4))/diameter;
                   % c_xydifp(counter)=distance(leftdisks(h,1:2)',centerpoint(1:2)')./diameter;
                    
                     c_xydifp(counter)=(leftdisks(h,1)-centerpoint(1)).^2;
                    c_xydifp(counter)=c_xydifp(counter)+(leftdisks(h,2)-centerpoint(2)).^2;
                   c_xydifp(counter)=(c_xydifp(counter)^.5)/diameter;

                    
                    newzdist(counter)=(j-pivot)*anisotropy/(diameter);
                    [m,minind]=min(distanceToPoint(rightdisks(:,1:2)',leftdisks(h,1:2)'));
                    goodmaximadifsp(counter)=(xydetdiameters(rightdisks(minind,4))-xydetdiameters(leftdisks(h,4)))./diameter;
                    goodmaximadifsp2(counter)=(xymaximavals(rightdisks(minind,4))-xymaximavals(leftdisks(h,4)))/intensity;
                    %xydisk_distancediff(counter)=(distance(rightdisks(minind,1:2)',leftdisks(h,1:2)'))./(diameter);
                    
                    xydisk_distancediff(counter)=(rightdisks(minind,1)-leftdisks(h,1)).^2;
                    xydisk_distancediff(counter)=xydisk_distancediff(counter)+(rightdisks(minind,2)-leftdisks(h,2)).^2;
                    xydisk_distancediff(counter)=(xydisk_distancediff(counter)^.5)/diameter;
       
                    counter=counter+1;
                end

            end
            if(j>pivot)
                for h=1:sizes2(1)
                    goodmaximasp(counter)=xymaximavals(rightdisks(h,4))./intensity;
                    gooddiamsp(counter)=xydetdiameters(rightdisks(h,4))/diameter;
                   % c_xydifp(counter)=distance(rightdisks(h,1:2)',centerpoint(1:2)')./diameter;
                    
                    
                    c_xydifp(counter)=(rightdisks(h,1)-centerpoint(1)).^2;
                    c_xydifp(counter)=c_xydifp(counter)+(rightdisks(h,2)-centerpoint(2)).^2;
                   c_xydifp(counter)=(c_xydifp(counter)^.5)/diameter;

                    
                    newzdist(counter)=(j+1-pivot)*anisotropy/(diameter);
                    [m,minind]=min(distanceToPoint(rightdisks(h,1:2)',leftdisks(:,1:2)'));
                    goodmaximadifsp(counter)=(xydetdiameters(rightdisks(h,4))-xydetdiameters(leftdisks(minind,4)))./diameter;
                    goodmaximadifsp2(counter)=(xymaximavals(rightdisks(h,4))-xymaximavals(leftdisks(minind,4)))/intensity;
                    
                    xydisk_distancediff(counter)=(rightdisks(h,1)-leftdisks(minind,1)).^2;
                    xydisk_distancediff(counter)=xydisk_distancediff(counter)+(rightdisks(h,2)-leftdisks(minind,2)).^2;
                    xydisk_distancediff(counter)=(xydisk_distancediff(counter)^.5)/diameter;
       
                    %xydisk_distancediff(counter)=(distance(rightdisks(h,1:2)',leftdisks(minind,1:2)'))./(diameter);
                    counter=counter+1;
                end
            end
        end
    end
end

%xy from center, diam, intensity, zdist, deltadiam deltaintensity diffxy
features=[c_xydifp;gooddiamsp;goodmaximasp;newzdist;goodmaximadifsp;goodmaximadifsp2;xydisk_distancediff];
