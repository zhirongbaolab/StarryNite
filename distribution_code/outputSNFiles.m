function outputSNfiles(directoryname,basename,timepoints,startime,endtime,downsample,ROIxmin,ROIymin)

%newDir=[directoryname,basename,'\matlabnuclei'];
newDir=[directoryname,basename,'_matlabnuclei'];
if(~(exist(newDir,'dir')))
    mkdir(newDir);
end


counter=1;
for time=startime:endtime

     outputfile=[newDir,'/matlabnuclei',num2str(time)];
     outputfile2=[newDir,'/matlabdiams',num2str(time)];
    fid=fopen(outputfile,'w');
        fid2=fopen(outputfile2,'w');
    if(~isempty(timepoints{counter})&&~isempty(timepoints{counter}.finalpoints))
        detpoints=round(timepoints{counter}.finalpoints);
       %detpoints=timepoints{counter}.finalaveragepoints;
        detpoints(:,1:2)=detpoints(:,1:2)./downsample;
        detpoints(:,1)=detpoints(:,1)+ROIxmin-1;
        detpoints(:,2)=detpoints(:,2)+ROIymin-1;
        count=fprintf(fid,'%d \t %d \t %d \n',detpoints');
     %  count=fprintf(fid,'%6.2f \t %6.2f \t %6.2f \n',detpoints');
        count=fprintf(fid2,'%d \n',timepoints{counter}.finaldiams'./downsample);
    end
    fclose(fid);
    fclose(fid2);
    counter=counter+1;

end %


