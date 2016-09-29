%given root filename including path, and timestep loads stack
function stack=loadNewCellStack(filename,slices, time,channel,ypixels)
i=1;
%slice=imread([filename,'-t',num2str(i,'%02d'),num2str(time,'%03d'),'-p','.tif']);
slice=imread([filename,'_s',num2str(i),'_t',num2str(time),'.tif']);
s=size(slice);
stack=ones(s(1),s(2),slices,'uint16');

stack(:,:,1)=slice;
for i=2:slices
    slice=imread([filename,'_s',num2str(i),'_t',num2str(time),'.tif']);
    stack(:,:,i)=slice;   
end
stack=reshape(stack,[size(slice),slices]);
if (channel==1)
    stack=stack(1:ypixels,1:ypixels,:);
else
    stack=stack(ypixels:ypixels*2,ypixels:ypixels*2,:);
end