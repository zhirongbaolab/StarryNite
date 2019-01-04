function convertToSideBySideImages(endtime,c1name,c2name,outname)
%note to simplify life by keeping all l/r images flipped imfipping these so
%standard acetree/sn paths/prameter files can be used

%revised newest nih for hobert
for time=1:endtime
        stack2=(loadSimpleStackTiff([c2name,'t',num2str(time),'.tif']));
        stack=(loadSimpleStackTiff([c1name,'t',num2str(time),'.tif']));
        s=size(stack);
        s(2)=s(2)*2;
        stackc=zeros(s,'uint16');
        stackc(:,1:s(2)/2,:)=fliplr(stack);
        stackc(:,s(2)/2+1:end,:)=fliplr(stack2);
    writeSimpleStackTiff([outname,'t',num2str(time),'.tif'],(stackc));
    
    %save([out_name,'_',num2str(c,'%04d'),'.mat'],'stack');
    
end

%{
endtime=99; 337;
c1name='L:\fanl\Imaging\For_Anthony\markers_from_yale\LF_DCR4315_ncs-1_BV293_08072017\LF_DCR4315_ncs-1_BV293_08072017_1_w2iSIM - TxRed - 600-50_s1_';
c2name='L:\fanl\Imaging\For_Anthony\markers_from_yale\LF_DCR4315_ncs-1_BV293_08072017\LF_DCR4315_ncs-1_BV293_08072017_1_w1iSIM - FITC - 525-50_s1_';
outname='L:\fanl\Imaging\For_Anthony\markers_from_yale\LF_DCR4315_ncs-1_BV293_08072017\LF_DCR4315_ncs-1_BV293_08072017_1_w2iSIM_combined_s1_';
convertToSideBySideImages(endtime,c1name,c2name,outname);


endtime= 337;
c1name='L:\fanl\Imaging\For_Anthony\markers_from_yale\LF_DCR6463_unc-86_1170_iSIM_11062017\LF_DCR6463_unc-86_1170_iSIM_11062017_1_w2iSIM - TxRed - 600-50_s3_';
c2name='L:\fanl\Imaging\For_Anthony\markers_from_yale\LF_DCR6463_unc-86_1170_iSIM_11062017\LF_DCR6463_unc-86_1170_iSIM_11062017_1_w1iSIM - FITC - 525-50_s3_';
outname='L:\fanl\Imaging\For_Anthony\markers_from_yale\LF_DCR6463_unc-86_1170_iSIM_11062017\LF_DCR6463_unc-86_1170_iSIM_11062017_1_w1iSIM_combined_s3_';
convertToSideBySideImages(endtime,c1name,c2name,outname);

for time=1:337
movefile([outname,'t',num2str(time),'.tif'],[outname,'_t',num2str(time),'.tif']);
end
%}
