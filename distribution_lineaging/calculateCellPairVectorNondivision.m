function [ datavector ] =calculateCellPairVectorNondivision(nucs,i,nucs2,j,anisotropy)
datavector=[    (nucs2.totalGFP(j)/nucs.totalGFP(i)),...
    nucs2.avgGFP(j)/nucs.avgGFP(i),...
distance(nucs.finalpoints(i,1:2)',nucs2.finalpoints(j,1:2)')/mean(nucs.selfdistance),...
distance(nucs.finalpoints(i,3)',nucs2.finalpoints(j,3)')*anisotropy(3)/mean(nucs.selfdistance)];

datavector(isinf(datavector))=0;
datavector(isnan(datavector))=0;
end

