%save greedy lineaging result in esequence data structure in  acetree readable format.
function saveGreedyNucleiFiles(esequence,time,directory,anisotropy,ROIxmin,ROIymin)
if(~exist('ROIxmin','var'))
    ROIxmin=1;
    ROIymin=1;
end
nextnum=1;%nuclear 'names'

%initialize traversed data structure
traversed=cell(time,1);
%for t=1:time
%    traversed{t}=zeros(size(esequence{t}.finalpoints,1),1,'uint8');
%end
maxpoint=false;
for t=1:time
    if(~isempty(esequence{t}.finalpoints))
        %note off by 2 correction here ROI is first included pixel so
        %pos=pos+roi-1  Acetree coordinate system is 0 origin rather than 1
        %which is the origin of second subtraction
        esequence{t}.finalpoints(:,1)=  esequence{t}.finalpoints(:,1)+ROIxmin-2;
        esequence{t}.finalpoints(:,2)=  esequence{t}.finalpoints(:,2)+ROIymin-2;
    end
end

%traverse points outputting info to array
%and then outputing to file once traversed

currentcounters=ones(1,time);
assembledData=cell(1,time);
traversed=cell(1,time);
for i=1:length(traversed)
    traversed{i}=zeros(size(esequence{i}.finalpoints,1),1,'uint8');
end


for t=1:time %we stop at time -1 because any initiation at endtime time is an isolated point
    for i=1:size(esequence{t}.finalpoints,1)
        
    %    if(esequence{t}.suc(i,1)~=-1&~esequence{t}.delete(i)) %not an isolated unlinked point or slated for deletion
       if(~esequence{t}.delete(i)) %not an isolated unlinked point or slated for deletion
    
            if(~traversed{t}(i)) %not traversed so its sucessors need to be entered
                %output the points info
                %not a death, or a dead successor or endtime
                if(t<time&esequence{t}.suc(i,1)~=-1&~esequence{esequence{t}.suc_time(i,1)}.delete(esequence{t}.suc(i,1)))
                    data=[-1,currentcounters(t+1)]; %pred,suc1
                    if(esequence{t}.suc(i,2)~=-1&~esequence{esequence{t}.suc_time(i,2)}.delete(esequence{t}.suc(i,2)))
                        data=[data,currentcounters(t+1)+1];%suc2
                    else
                        data=[data,-1];%suc2
                    end
                else
                    data=[-1,-1,-1]; %pred,suc1
                end
                
                
                data=[data,esequence{t}.finalpoints(i,:),esequence{t}.finaldiams(i),nextnum,esequence{t}.totalGFP(i)];%x,y,z,diam,nextnum(used in name)
                
                %assembledData{t}(currentcounters(t))=data;
                assembledData{t}=[assembledData{t};data];
                
                
                currentcounters(t)=currentcounters(t)+1;
                traversed{t}(i)=1;
                nextnum=nextnum+1;
                if(t<time&esequence{t}.suc(i,1)~=-1&~esequence{esequence{t}.suc_time(i,1)}.delete(esequence{t}.suc(i,1)))
                    %recurse daughter 1
                    [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,t,i,1,nextnum,time);
                    if(esequence{t}.suc(i,2)~=-1&~esequence{esequence{t}.suc_time(i,2)}.delete(esequence{t}.suc(i,2)))
                        [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,t,i,2,nextnum,time);
                    end
                end
            end
        end
    end
end




%ouput assembled data

for i=1:time%end_time %per time
    
    outputfile=[directory,'/t',num2str(i,'%03d'),'-nuclei'];
    fid=fopen(outputfile,'w');
    
    
    for p=1:size(assembledData{i},1)%per point
  %note arbitrary attempt to keep gfp total in grange acetree will read
            
        fprintf(fid,'%d,  %d,  %d,  %d,  %d,  %d, %d,  %2.1f,  %d,  %s,  %d, %s\n',...
            p,1,round(assembledData{i}(p,1:7)),['Nuc',num2str(assembledData{i}(p,8))],uint16(assembledData{i}(p,9)/256),' 0, 0, 0, , 0, 0, 0, 0, 0,');

    end
    
    
    fclose(fid);

end

end %end function

function [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,t,i,daughternum,nextnum,time)
%fill this in next above should be good
    daughteri=esequence{t}.suc(i,daughternum);
    daughteri_t=esequence{t}.suc_time(i,daughternum);
    
    %insert here interpolation output
    %not exectuted in single frame jump
    tlocal=t+1;
    while(tlocal<daughteri_t&tlocal<=time)
        timespan=daughteri_t-t;
        pointstart=esequence{t}.finalpoints(i,:);
        pointend=esequence{daughteri_t}.finalpoints(daughteri,:);
        currentpoint=pointstart*(1-(tlocal-t)/timespan)+pointend*(tlocal-t)/timespan;
        %if interpolating into endtime (only an event that happens in
        %training) dont put a successor in endtime
        if(tlocal<time)
            data=[currentcounters(tlocal-1)-1,currentcounters(tlocal+1),-1]; %pred,suc1,suc2
        else
             data=[currentcounters(tlocal-1)-1,-1,-1]; %pred,suc1,suc2
        end
       % data=[data,(esequence{t}.finalpoints(i,:)),esequence{t}.finaldiams(i),nextnum,esequence{t}.totalGFP(i)];%x,y,z,diam,nextnum(used in name)
        data=[data,(currentpoint),esequence{t}.finaldiams(i),nextnum,esequence{t}.totalGFP(i)];%x,y,z,diam,nextnum(used in name)

        assembledData{tlocal}=[assembledData{tlocal};data];
        currentcounters(tlocal)=currentcounters(tlocal)+1;
        tlocal=tlocal+1;
       %nothing traversed during interpolation traversed{daughteri_t}(daughteri)=1;
        nextnum=nextnum+1;
    end
    %if a FN has walked us over endtime dont output and done
    if(daughteri_t<=time)
        if(daughteri_t+1<=time)
            daughterisuc1=esequence{daughteri_t}.suc(daughteri,1);
            %output the actual points info
            if(daughterisuc1~=-1)
                data=[currentcounters(daughteri_t-1)-1,currentcounters(daughteri_t+1)]; %pred,suc1
                if(esequence{daughteri_t}.suc(daughteri,2)~=-1)
                    data=[data,currentcounters(daughteri_t+1)+1];%suc2
                else
                    data=[data,-1];%suc2
                end
            else
                data=[currentcounters(daughteri_t-1)-1,-1,-1]; %pred,suc1 %changed this from currentcounters(t)-1
            end
        else
            data=[currentcounters(daughteri_t-1)-1,-1,-1];
        end
        data=[data,esequence{daughteri_t}.finalpoints(daughteri,:),esequence{daughteri_t}.finaldiams(daughteri),nextnum,esequence{daughteri_t}.totalGFP(daughteri)];%x,y,z,diam,nextnum(used in name)
        
        
        assembledData{daughteri_t}=[assembledData{daughteri_t};data];;
        currentcounters(daughteri_t)=currentcounters(daughteri_t)+1;
        traversed{daughteri_t}(daughteri)=1;
        nextnum=nextnum+1;
        %having output this cell output its daughters if exist
        if(esequence{daughteri_t}.suc(daughteri,1)~=-1&&daughteri_t<time&~esequence{esequence{daughteri_t}.suc_time(daughteri,1)}.delete(esequence{daughteri_t}.suc(daughteri,1)))
            [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,daughteri_t,daughteri,1,nextnum,time);
        end
        if(esequence{daughteri_t}.suc(daughteri,2)~=-1&&daughteri_t<time&~esequence{esequence{daughteri_t}.suc_time(daughteri,2)}.delete(esequence{daughteri_t}.suc(daughteri,2)))
            [traversed,currentcounters,assembledData,nextnum]=recursive_traverse_output(esequence,traversed,currentcounters,assembledData,daughteri_t,daughteri,2,nextnum,time);
        end
    end
end
