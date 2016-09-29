function [merges,ranges]=filter_boundary(centers,xyzmax,ranges,logodds,overlaps,xymax,xymaximavals,xydetdiameters,anisotropy,xycoverage)
%test every pair of nuclei that overlap by a disk
%return list of scoreing of iv they need to be merged
%and updated ranges for
merges=[];
rangesnew=ranges;
for i=1:length(overlaps)
    ol=overlaps{i};
    if (~isempty(ol))%overlap exists
        for h=1:length(ol)
            [s1,s2,rangesnew{i},rangesnew{ol(h)}]=mergeDecision(i,ol(h),xyzmax(i,:),xyzmax(ol(h),:),centers,ranges,logodds,xymax,xymaximavals,xydetdiameters,anisotropy,xycoverage);
            merges=[merges;i,ol(h),s1,s2];
        end
    end
end
ranges=rangesnew;