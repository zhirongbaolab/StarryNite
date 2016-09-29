function filtered=tiledogfilter(data,sigma,anisotropy)


%chunkxy=150;%size of chunk we can handle in memory
%chunkz=21;%21 is problem


chunkxy=500;%size of chunk we can handle in memory
chunkz=100;%21 is problem

maxpaddingxy=floor(sigma*1.6/2); %ammount of tile that will be invalid bc filter will run off edge
maxpaddingz=floor(sigma*1.6/anisotropy/2);

s=size(data);

filtered=zeros(s,'single');
for cx=1:chunkxy:s(1)
    for cy=1:chunkxy:s(2)
        for cz=1:chunkz:s(3)
            %padded out so whole chunk is good
            dchunk=single(data(max(1,cx-maxpaddingxy):min(s(1),cx+chunkxy+maxpaddingxy),max(1,cy-maxpaddingxy):min(s(2),cy+chunkxy+maxpaddingxy),max(1,cz-maxpaddingz):min(s(3),cz+chunkz+maxpaddingz)));
            %filter 
            %FFT dogfilter now uses replication padding, so edge wrapping
            %should not be a problem anymore
            dchunk=dogfilter(dchunk,sigma,anisotropy);
            
            schunk=size(dchunk);
            %cut away contaminated edge
            if(cx==1) 
                goodstartx=1;
            else
                goodstartx=min(maxpaddingxy,cx)+1;
            end
            if(cy==1) 
                goodstarty=1;
            else
                goodstarty=min(maxpaddingxy,cy)+1;
            end
            if(cz==1 ) 
                goodstartz=1;
            else
                %if padding is bigger than increment need to 
                goodstartz=min(maxpaddingz,cz)+1;
            end
            
            if(cx+chunkxy>=s(1)) 
                goodendx=schunk(1);
            else
                goodendx=schunk(1)-maxpaddingxy;
            end
            if(cy+chunkxy>=s(2)) 
                goodendy=schunk(2);
            else
                goodendy=schunk(2)-maxpaddingxy;
            end
            if(cz+chunkz>=s(3)) 
                goodendz=schunk(3);
            else
               
                goodendz=schunk(3)-min(maxpaddingz,s(3)-(cz+chunkz)+1);
            end
            dchunk=dchunk(goodstartx:goodendx,goodstarty:goodendy,goodstartz:goodendz);
            schunk=size(dchunk);
            filtered(cx:cx+schunk(1)-1,cy:cy+schunk(2)-1,cz:cz+schunk(3)-1)=dchunk;
        end %dim 1 loop
    end  %dim 2 loop
end%dim 3 loop