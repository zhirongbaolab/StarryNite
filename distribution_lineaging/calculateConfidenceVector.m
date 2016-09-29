function [ confidencedata ] = calculateConfidenceVector( nuclei,parameters )
%calculates matrix of confidence cue values given timepoint structure
%%sumgfp, nndistance (diams) z, xy, aspectratio, logodds_avg, avg gfp

confidencedata=zeros(size(nuclei.finalpoints,1),6);
[integratedGFP,area]=integrateGFP(nuclei,parameters);
%confidencedata(:,1)=log(integratedGFP+1);
confidencedata(:,1)=integratedGFP./mean(integratedGFP);
distances=distance_anisotropic(nuclei.finalpoints',nuclei.finalpoints',parameters.anisotropyvector);

distancesz=distance(nuclei.finalpoints(:,3)',nuclei.finalpoints(:,3)').*parameters.anisotropyvector(3);
distancesxy=distance(nuclei.finalpoints(:,1:2)',nuclei.finalpoints(:,1:2)');

if(size(nuclei.finalpoints,1)>1)
    for i=1:size(nuclei.finalpoints,1)
        distances(i,i)=Inf;
        %distancesxy(i,i)=Inf;
    end
end
[m,ivals]=min(distances);

confidencedata(:,4)=nuclei.aspectratio./mean(nuclei.aspectratio);
%this is a work around for a bug in earlier versions of detection algorithm
%allows old .mat files to be relineaged
confidencedata(:,5)=nuclei.mergedlogoddssum(1:size(nuclei.finalpoints,1));
confidencedata(:,6)=integratedGFP./area;
confidencedata(area==0,6)=0;
if(mean(confidencedata(:,6))~=0)
    confidencedata(:,6)=confidencedata(:,6)./mean(confidencedata(:,6));
end
for i=1:size(nuclei.finalpoints,1)
    confidencedata(i,5)= confidencedata(i,5)./length(nuclei.merged_sliceindicies{i});
    
    confidencedata(i,2)=log((distancesz(i,ivals(i)))./(mean(nuclei.selfdistance))+1);
    confidencedata(i,3)=log((distancesxy(i,ivals(i)))./(mean(nuclei.selfdistance))+1);
    
end
if(mean(confidencedata(:,5))~=0)
    confidencedata(:,5)=confidencedata(:,5)./mean(confidencedata(:,5));
end

end

