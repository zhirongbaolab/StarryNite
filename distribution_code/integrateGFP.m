function [integratedGFP,area,maxslice,diskmax]=integrateGFP(dataset,parameters)
integratedGFP=zeros(length(dataset.merged_sliceindicies),1);
area=zeros(length(dataset.merged_sliceindicies),1);
maxslice=zeros(length(dataset.merged_sliceindicies),1);
diskmax=zeros(length(dataset.merged_sliceindicies),1);
 for n=1:length(dataset.merged_sliceindicies)  
   
            slices=dataset.merged_sliceindicies{n};
%            filterslice=true;
            
 %           if(filterslice)
            
                maxvalue=max(dataset.diskintensity(slices));
                slices=slices(dataset.diskintensity(slices)>=maxvalue*parameters.boundary_percent);
                maxslice(n)=maxvalue;
                %          end
           % diskmax(n)=prctile(dataset.diskMax(slices),75);
            for s=1:length(slices)
                 area(n)=area(n)+dataset.diskArea(slices(s));
                integratedGFP(n)=integratedGFP(n)+dataset.diskGFPsums(slices(s));
               diskmax(n)=max(diskmax(n),dataset.diskMax(slices(s)));
            end  

end
   
       

        