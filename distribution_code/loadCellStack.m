%given root filename including path, and timestep loads stack
function stack=loadCellStack(filename,slices, time)
i=1;
slice=imread([filename,'-t',num2str(time,'%03d'),'-p',num2str(i,'%02d'),'.tif']);
s=size(slice);
stack=ones(s(1),s(2),slices,'uint8');

stack(:,:,1)=slice;
for i=2:slices
    slice=imread([filename,'-t',num2str(time,'%03d'),'-p',num2str(i,'%02d'),'.tif']);
    stack(:,:,i)=slice;   
end
stack=reshape(stack,[size(slice),slices]);
