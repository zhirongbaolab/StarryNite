function [integratedGFP,area]=integrateGFP(dataset,parameters)
integratedGFP=zeros(length(dataset.merged_sliceindicies),1);
area=zeros(length(dataset.merged_sliceindicies),1);

 for n=1:length(dataset.merged_sliceindicies)  
   
            slices=dataset.merged_sliceindicies{n};
            filterslice=true;
            
            if(filterslice)
                maxvalue=max(dataset.diskintensity(slices));
                slices=slices(dataset.diskintensity(slices)>=maxvalue*parameters.boundary_percent);
            end
      
            for s=1:length(slices)
                 area(n)=area(n)+dataset.diskArea(slices(s));
                integratedGFP(n)=integratedGFP(n)+dataset.diskGFPsums(slices(s));
            end  

end
   
       

        