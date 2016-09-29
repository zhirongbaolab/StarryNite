function  outputAceTreeSlice(X,directoryname, basename,time,minvalue, maxvalue, downsample,red)
%newDir=[directoryname,basename,'\tif'];
%X=im2uint8(X);
if (red)
    newDir=[directoryname,'image/tifR'];
else
newDir=[directoryname,'image/tif'];
end
if(~(exist(newDir,'dir')))
    mkdir(newDir);
end

s=size(X);
%minx=min(min(min(X)));
%maxx=max(max(max(X)));
%X=single(X);
%X=(X-minx)/(maxx-minx);
%X=(X-single(minvalue))./single(maxvalue-minvalue);
%maxthere=max(max(max(X)));
for slice=1:s(3)
   % if(s(3)<100)
   s_out=X(:,:,slice);
   %s_out=single(s_out);
   s_out=((single(s_out)-single(minvalue))./single(maxvalue-minvalue));
   s_out(s_out<0)=0;
   s_out(s_out>1)=1;
   s_out=uint8(255*s_out);
        imwrite(s_out,sprintf('%s/%s-t%03d-p%02d.tif',newDir,basename,time,slice),'tiff','Compression','none');
    %else
    %    imwrite(((X(:,:,slice))),sprintf('%s/%s-t%03d-p%03d.tif',newDir,basename,time,slice),'tiff','Compression','none');
    %end
end
