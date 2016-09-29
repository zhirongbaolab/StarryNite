function [ value] = approximateMatrixPercentile(data,percentile,resolution)
% compute approximate percentile of matrix data using histogram of resolution resolution
%bin values and interpolate within bin to find rough percentile 
minval=min(min(min(data)));
maxval=max(max(max(data)));

bins=linspace(double(minval),double(maxval),resolution);

histo=histc(reshape(data(:,:,1),1,size(data,1)*size(data,2)),bins);
for z=2:size(data,3)
    histo=histo+histc(reshape(data(:,:,z),1,size(data,1)*size(data,2)),bins);
end
%index of item corresponding to percentile
desired=percentile/100*prod(size(data));
%now walk along histogram
below=0;
ind=1;
while below<desired
    below=below+histo(ind);
    if(below<desired)
        ind=ind+1;
    end
end
%ind contains index
 prior=sum(histo(1:ind-1));
 additional=desired-prior;
 weight=additional/histo(ind);
 if(ind==1) 
     value=bins(ind)*(1-weight)+bins(ind+1)*weight;
 else
      value=bins(ind-1)*(1-weight)+bins(ind)*weight;
 end

end

