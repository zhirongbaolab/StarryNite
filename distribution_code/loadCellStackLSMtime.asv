function stack=loadCellStackLSMtime(name,time,channel,nPlanes)
%hopefully general lsm reader
%loads stack that is particular timepoint and channel from larger lsm file
% size of stack is min of nPlanes and number of planes from header
% channel is overridden to chanel 1 if only one channel
 lsminformation=lsminfo(name);
 
W=lsminformation.DimensionX;
L=lsminformation.DimensionY;
nPlanes=min(nPlanes,lsminformation.DimensionZ);
stack=zeros(L,W,nPlanes);
channel=min(channel,lsminformation.NUMBER_OF_CHANNELS);

for j=1:nPlanes

    %embedded in this is the assumption which has been case for every lsm
    %file I've seen that it interleaves them main file and thumbnails
    subimagenum=((time-1)*2)*nPlanes+(j-1)*2+1;
    %im=imread(name,'Index',subimagenum);
    %for some reason imread no longer works
    %revert to non matlab reader
    im=tiffread32(name,subimagenum);
    im=im.data;
    stack(:,:,j)=im(:,:,channel);

end

