%parameterfile assumed to be defined and point to parameter file  

%paramater file is formated as a series of matlab commands

fid = fopen(parameterfile, 'rt');

while feof(fid) == 0
   tline = fgetl(fid);
   if(tline~=-1) %returns -1 if trailing hard returns
        eval(tline);
   end
end
fclose(fid);