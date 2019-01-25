%save without lineaging in acetree readable format.
%

%base='l:/santella/mouse/Anthony/paperdata/';
%base='l:/santella/zebrafish/paperdata/';
%base='l:/santella/keller_d/paperdata/';
%base='l:/santella/NIH_data/Embryo_20101001/deconresults/annot/dats/nuclei/';
counter=1;
for i=start_time:end_time
   % e=esequence{i-start_time+1};
   e=esequence{i}; %I'm not sure why the above (whcih doesnt work) was set if I changed this behavior at some point and did so inconsistently 4/29/2015 -as
    outputfile=[base,'t',num2str(i,'%03d'),'-nuclei'];
    fid=fopen(outputfile,'w');
    

        %{   
            bounds=[90,610,400,580,35,60];%fly expanded2
            %stack=X(bounds(1):bounds(2),bounds(3):bounds(4),bounds(5):bounds(6));
            
       %   bounds=[0,460,170,260,25,45];% a mouse expanded
 %bounds=[1,460,170,260,25,45];% a mouse expanded
       
               bounds=[60,202,-1,202,0,12];%zebrafish late 
               bounds=[-1,600,-1,600,0,25];%zebrafish early no bounds
               
                 bounds=[90,610,400,590,35,60];%fly expanded3
                  bounds=[90,610,420,560,35,60];%fly expanded
%}    
    if(ROI)
       % note off by 2 correction here ROI is first included pixel so
        %pos=pos+roi-1;  Acetree coordinate system is 0 origin rather than 1
        %which is the origin of second subtraction
        if( ~isempty(e.finalpoints))
        e.finalpoints(:,1)=e.finalpoints(:,1)+ROIxmin-2;
        e.finalpoints(:,2)=e.finalpoints(:,2)+ROIymin-2;
        end
    else
           if( ~isempty(e.finalpoints))  
         e.finalpoints(:,1)=e.finalpoints(:,1)-1;
        e.finalpoints(:,2)=e.finalpoints(:,2)-1;
           end
    end
outcounter=1;
    for j=1:length(e.finaldiams)
    
        pcurrent=e.finalpoints(j,:);
        %test bounds bc acetree editing seems broken on large data sets
       % if(pcurrent(1)>bounds(1)&pcurrent(1)<bounds(2)&pcurrent(2)>bounds(3)&pcurrent(2)<bounds(4)&pcurrent(3)>bounds(5)&pcurrent(3)<bounds(6))
     
            fprintf(fid,'%d, \t %d, \t %d, \t %d, \t %d, \t %d, %d, \t %d, \t %d, \t %s, \t %d, \n',outcounter,1,-1,-1,-1,int16(e.finalpoints(j,1)),int16(e.finalpoints(j,2)),int16(e.finalpoints(j,3)),e.finaldiams(j),['Nuc',num2str(j*i)],0);
            outcounter=outcounter+1;
      %  end
    end
    fclose(fid);
    counter=counter+1;
end