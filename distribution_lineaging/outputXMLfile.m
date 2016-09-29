function     outputXMLfile(xmlname,xyres,zres,slices,zipname,embryodir,embryonumber,newimage);

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

file=fopen(xmlname,'w');
fprintf (file, '<?xml version=''1.0'' encoding=\''utf-8\''?>\n');
fprintf (file, '<embryo>\n');
if (newimage)
    %if looking at new images delete old images after use
   % rmdir([embryodir,'image/'],'s');
    fprintf (file,'<useStack type="1"/>\n');
    fprintf (file,['<image file="',embryodir,embryonumber,'_t1.tif"/>\n']);
else
    fprintf (file,['<image file="',embryodir,'image/tif/',embryonumber,'-t001-p01.tif"/>\n']);
end
fprintf (file,['<nuclei file="',zipname,'"/>\n']);

fprintf (file,'<end index="475"/>\n');
fprintf(file,['<resolution xyRes="',num2str(xyres),'" zRes="',num2str(zres),'" planeEnd="',num2str(slices),'"/> <exprCorr type="blot"/>\n']); 
fprintf (file,'<polar size="15"/>\n');
fprintf (file,'</embryo>');
fclose(file);
end

