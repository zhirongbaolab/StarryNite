function stack=loadCellStackLSM(name,time,channel,nPlanes)
%name is filename of combined tif
%channel is which r/g to get out 2 is nuclei green 1 is red
%not botherting to write a general loader assuming all the mouse images are
%paired brightfield? and 2 channel florescence images

block=1;
timepoint=1;

for t=1:time
    filename=[name,'_GR',num2str(timepoint),'_B',num2str(block) ,'.lsm'];
    if (~exist(filename,'file'))
        block=block+1;
        timepoint=1;
         filename=[name,'_GR',num2str(timepoint),'_B',num2str(block) ,'.lsm'];
         if (~exist(filename,'file'))
             error(['filename ',filename,' not found']);
         end
    else
        if(t~=time)
            timepoint=timepoint+1;
        end
    end
end

filename=[name,'_GR',num2str(timepoint),'_B',num2str(block) ,'.lsm'];

% if exists file then imrement till time if not 
%block is final block #
%timepoint is timepoint in block

imageinfo=imfinfo(filename);
L=imageinfo(1).Height;

W=imageinfo(1).Width;
depth=length(imageinfo)/2;

stack=zeros(L,W,min(nPlanes,depth),'uint8');
for j=1:depth
    subimagenum=(j-1)*2+1;
    %subimagenum=((time-1)*2)*nPlanes+(j-1)*2+1;
    im=imread(filename,subimagenum);
    stack(:,:,j)=im(:,:,1);

end

%stack=stack./(2^16);