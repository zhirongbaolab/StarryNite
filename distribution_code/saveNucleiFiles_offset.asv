%save without lineaging in acetree readable format.
function saveNucleiFiles_offset(tracksFinal,esequence,time,directory,anisotropy,abs_start_time)

%preprocessing hack if nan fill in with previous result (ie stands still
%during gaps instead of interpolating


%leave past end nan
for ct=1:size(tracksFinal,1)
       for st=1:size(tracksFinal(ct).tracksCoordAmpCG,1)
           started=false;
           endtime=max(find(~isnan(tracksFinal(ct).tracksCoordAmpCG(st,:))));
            for p=1:8:size(tracksFinal(ct).tracksCoordAmpCG,2)
                if(~isnan(tracksFinal(ct).tracksCoordAmpCG(st,p)))
                    started=true;
                end
                if(started&&p<=endtime&& isnan(tracksFinal(ct).tracksCoordAmpCG(st,p)))
                    tracksFinal(ct).tracksCoordAmpCG(st,p:p+3)=tracksFinal(ct).tracksCoordAmpCG(st,p-8:p-5);
                end
            end
       end
end

    
    
currentcounters=ones(1,time);
assembledData=cell(1,time);

 abs_subtrack_number=1;
for ct=1:size(tracksFinal,1)%iterate over all compound tracks
    %absoute start time of this compound track
    cts_start=tracksFinal(ct).seqOfEvents(1,1);
    
    %initialize an array of pointers to final resting place of each nuc in
    %output array
    tracksFinal(ct).RowNumbers=zeros(size(tracksFinal(ct).tracksFeatIndxCG));
    
    for st=1:size(tracksFinal(ct).seqOfEvents,1) %iterate over events in compound tracks
     
        current_subtrack=tracksFinal(ct).seqOfEvents(st,3);
        
        if(tracksFinal(ct).seqOfEvents(st,2)==1)% track born add it to assembled data
            % as long as track exists output it
            localtp=1;
            abs_starttime=tracksFinal(ct).seqOfEvents(st,1);
           % local_endtime=max(find(~isnan(tracksFinal(ct).tracksCoordAmpCG(st,:))));
            
            local_index=(abs_starttime-cts_start+localtp-1) *8+1;
            abs_currenttime=abs_starttime+localtp-1;
            local_index_single=abs_starttime-cts_start+localtp-1;
            
            %predecessor/suc for division
                
            %while local_index<local_endtime&&
            while(local_index_single<size(tracksFinal(ct).tracksFeatIndxCG,2)&&~isnan(tracksFinal(ct).tracksCoordAmpCG(current_subtrack,local_index)))
                
                assembledData{abs_currenttime}(currentcounters(abs_currenttime),1:4)=tracksFinal(ct).tracksCoordAmpCG(current_subtrack,local_index:local_index+3);
                 assembledData{abs_currenttime}(currentcounters(abs_currenttime),5)=abs_subtrack_number;
                
                 %put in diam
                if(tracksFinal(ct).tracksFeatIndxCG(current_subtrack,local_index_single+1)~=0)
                 assembledData{abs_currenttime}(currentcounters(abs_currenttime),6)=esequence{abs_currenttime+abs_start_time}.finaldiams(tracksFinal(ct).tracksFeatIndxCG(current_subtrack,local_index_single+1)); 
                else
                 assembledData{abs_currenttime}(currentcounters(abs_currenttime),6)= assembledData{abs_currenttime-1}(tracksFinal(ct).RowNumbers(current_subtrack,local_index_single),6);
                end
                
                 assembledData{abs_currenttime}(currentcounters(abs_currenttime),7:9)=[-1,-1,-1];
                  
                 %predecessor
                if(localtp>1)
                    assembledData{abs_currenttime}(currentcounters(abs_currenttime),7)=tracksFinal(ct).RowNumbers(current_subtrack,local_index_single);
                    assembledData{abs_currenttime-1}(tracksFinal(ct).RowNumbers(current_subtrack,local_index_single),8)=currentcounters(abs_currenttime);
                 
                else
                    if(~isnan( tracksFinal(ct).seqOfEvents(st,4)))
                         assembledData{abs_currenttime}(currentcounters(abs_currenttime),7)=tracksFinal(ct).RowNumbers(tracksFinal(ct).seqOfEvents(st,4),local_index_single);
                         assembledData{abs_currenttime-1}(tracksFinal(ct).RowNumbers(tracksFinal(ct).seqOfEvents(st,4),local_index_single),9)=currentcounters(abs_currenttime);
                     end
                end
                 
             
                 
                %save where put the current cell for future reference when
                %divisions
                tracksFinal(ct).RowNumbers(current_subtrack,local_index_single+1)=currentcounters(abs_currenttime);
                
                currentcounters(abs_currenttime)=currentcounters(abs_currenttime)+1;%inc next open space
                localtp=localtp+1;
                local_index=(abs_starttime-cts_start+localtp-1) *8+1;
                abs_currenttime=abs_starttime+localtp-1;
                local_index_single=abs_starttime-cts_start+localtp-1;
            
            end
            abs_subtrack_number= abs_subtrack_number+1;
            
        end
    end
    
end
'ok'


end_time=length(tracksFinal(1).tracksCoordAmpCG)/8;

for i=1:time%end_time %per time
    
    outputfile=[directory,'/t',num2str(i+abs_start_time,'%03d'),'-nuclei'];
    fid=fopen(outputfile,'w');
    
    
    for p=1:size(assembledData{i},1)%per point
        
        %for moment divisions are not output******
        fprintf(fid,'%d, \t %d, \t %d, \t %d, \t %d, \t %0d, %0d, \t %3.3f, \t %d, \t %s, \t %d, \n',p,1,assembledData{i}(p,7),assembledData{i}(p,8),assembledData{i}(p,9),int16(assembledData{i}(p,1)),int16(assembledData{i}(p,2)),assembledData{i}(p,3)./anisotropy,assembledData{i}(p,6),['Nuc',num2str(assembledData{i}(p,5))],uint16(assembledData{i}(p,4)));

    end
    
    
    fclose(fid);

end

%{
%totally screwed up

counter=1;

end_time=length(tracksFinal(1).tracksCoordAmpCG)/8;
for i=1:end_time %per time
    
    outputfile=[directory,'/t',num2str(i,'%03d'),'-nuclei'];
    fid=fopen(outputfile,'w');
    outcounter=1;
    
    for t=1:length(tracksFinal) %per track
        trackarray=tracksFinal(t).tracksCoordAmpCG;
        
        for p=1:size(trackarray,1)%per point within track
            diameter=10;
            if(~isnan(trackarray(p,1+8*(t-1))))
                %for moment divisions are not output******
                fprintf(fid,'%d, \t %d, \t %d, \t %d, \t %d, \t %d, %d, \t %d, \t %d, \t %s, \t %d, \n',outcounter,1,outcounter,outcounter,-1,trackarray(p,1+8*(t-1)+1),trackarray(p,1+8*(t-1)+2),trackarray(p,1+8*(t-1)+3),diameter,['Nuc',num2str(counter)],trackarray(p,1+8*(t-1)+4));
                outcounter=outcounter+1;
            else%invalid nuc hack for moment, ** just causes inefficiency on disk
                fprintf(fid,'%d, \t %d, \t %d, \t %d, \t %d, \t %d, %d, \t %d, \t %d, \t %s, \t %d, \n',outcounter,-1,-1,-1,-1,0,0,0,0,['Nuc',num2str(counter)],0);
                outcounter=outcounter+1;
            end
        end
        
        
    end
    fclose(fid);
    counter=counter+1;
end
%}