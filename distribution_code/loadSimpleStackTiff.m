function stack=loadSimpleStackTiff(name)



imageinfo=imfinfo(name);
L=imageinfo(1).Height;

W=imageinfo(1).Width;
depth=length(imageinfo);

stack=zeros(L,W,depth);
for j=1:depth
    im=imread(name,j);
    stack(:,:,j)=im;

end

