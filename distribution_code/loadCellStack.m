%given root filename including path, and timestep loads stack
function stack=loadCellStack(filename,slices, time,delim1,delim2)
if(~exist('delim1'))
    delim1='-';
end

if(~exist('delim2'))
    delim2='-';
end

i=1;
%slice=imread([filename,delim1,'t',num2str(time,'%03d'),delim2,'p',num2str(i,'%02d'),'.tif']);
slice=imread([filename,delim1,'t00',num2str(time),delim2,'p',num2str(i,'%02d'),'.tif']);

s=size(slice);
%stack=ones(s(1),s(2),slices,class(slice));
'hack fix operational'
stack=ones(260,260,slices,class(slice));

stack(1:s(1),1:s(2),1)=slice;
for i=2:slices
%    slice=imread([filename,delim1,'t',num2str(time,'%03d'),delim2,'p',num2str(i,'%02d'),'.tif']);
     slice=imread([filename,delim1,'t00',num2str(time),delim2,'p',num2str(i,'%02d'),'.tif']);
    [s1,s2]=size(slice);
    stack(1:s1,1:s2,i)=slice;   
end
stack=reshape(stack,[size(slice),slices]);
