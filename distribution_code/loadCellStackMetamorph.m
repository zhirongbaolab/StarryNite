function stack=loadCellStackMetamorph(basename,time, channel,nPlanes,offset,zeropadding)
%name is base filename of combined tif
%time is timepoint to load
%channel is which r/g to get out 2 is nuclei green 1 is red
%L is image length
%totalW is side by side slice image width
%offset is an [x,y] vector of integers representing the number of pixels to offset images by
%positive x shifts right image further to the right
%positive y shifrs right image down
%embryoIm is a 2x2 matrix specifing the part of the image containing the
%embryo.  It is of the form, [first vertical pixel, last vertical pixel;
%first horizontal pixel, last horizontal pixel];

%{
name='L:\new embryos until twitching\original metamorph data\s2\BV24_wt_061710_s2_t001.TIF';
%test=imread(name);

test=loadCellStackMetamorph(name,2,512,1024,30,[0,0]);
%}
%L = 512;
%totalW = 1024;

%nPlanes = 30;

%parameters for images size, may need to be changed

if (zeropadding)
    name=[basename,'_t',num2str(time,'%03d'),'.TIF'];
else
    name=[basename,'_t',num2str(time),'.TIF'];
end

imageinfo=imfinfo(name);
nPlanes=min(nPlanes,length(imageinfo));%override specified number of slices if it is more than those present
L=imageinfo(1).Height;

totalW=imageinfo(1).Width;


W = totalW/2;
%set parts of sub image you will copy taking offsets into account
subIm=cell(1,2);
subIm{1}=cell(1,2);
subIm{2}=cell(1,2);
%if nargin>=7
%    subIm{1}={[embryoIm(1,1) embryoIm(1,2)],...
%        [embryoIm(2,1) embryoIm(2,2)]};
%else
    subIm{1}={[1,L],[1,W]};
%end
Lstart=zeros(1,2);Lfinish=zeros(1,2);
Wstart=zeros(1,2);Wfinish=zeros(1,2);
Lstart(1)=1;Lfinish(1)=subIm{1}{1}(2)-subIm{1}{1}(1)+1;
Wstart(1)=1;Wfinish(1)=subIm{1}{2}(2)-subIm{1}{2}(1)+1;


%shifts right image relative to left image
if offset(2)>=0
    Lstart(2)=1+offset(2);
    Lfinish(2)=subIm{1}{1}(2)-subIm{1}{1}(1)+1;
    if offset(1)>=0
        subIm{2}={[subIm{1}{1}(1),subIm{1}{1}(2)-offset(2)],...
            [subIm{1}{2}(1)+W,W+subIm{1}{2}(2)-offset(1)]};
        Wstart(2)=1+offset(1);
        Wfinish(2)=subIm{1}{2}(2)-subIm{1}{2}(1)+1;
    else
        subIm{2}={[subIm{1}{1}(1),subIm{1}{1}(2)-offset(2)],...
            [subIm{1}{2}(1)+W-offset(1),W+subIm{1}{2}(2)]};
        Wstart(2)=1;
        Wfinish(2)=abs(subIm{1}{2}(2)-subIm{1}{2}(1))+offset(1)+1;
    end        
else
    Lstart(2)=1;
    Lfinish(2)=subIm{1}{1}(2)-subIm{1}{1}(1)+offset(2)+1;
    if offset(1)>=0
        subIm{2}={[subIm{1}{1}(1)-offset(2),subIm{1}{1}(2)],...
            [subIm{1}{2}(1)+W,W+subIm{1}{2}(2)-offset(1)]};
        Wstart(2)=1+offset(1);
        Wfinish(2)=subIm{1}{2}(2)-subIm{1}{2}(1)+1;
    else
        subIm{2}={[subIm{1}{1}(1)-offset(2),subIm{1}{1}(2)],...
            [subIm{1}{2}(1)+W-offset(1),W+subIm{1}{2}(2)]};
        Wstart(2)=1;
        Wfinish(2)=subIm{1}{2}(2)-subIm{1}{2}(1)+offset(1)+1;
    end  
end

stack=zeros(L,W,nPlanes);
for j=1:nPlanes
    k=channel;
        if (nargin>=7)
            stack(Lstart(k):Lfinish(k),Wstart(k):Wfinish(k),j)=...
                imread(name,'tiff',j,'PixelRegion',subIm{k});
        else
            stack(Lstart(k):Lfinish(k),Wstart(k):Wfinish(k),j)=...
                imread(name,'tiff',j,'PixelRegion',subIm{k});
        end
end

%stack=stack./(2^16);