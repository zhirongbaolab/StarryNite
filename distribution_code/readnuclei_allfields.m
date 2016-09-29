function [M cellNames forcedCellNames] = readnuclei_allfields(filename, col)
%CREATE DATE: apr 21, 2008
%Reads a nuclei file 
%with the following structure :
%1)file specific index (int)
%2)is valid (int)
%3)predessesor (int)
%4)successor1 or -1 (int)
%5)successor2 or -1 (int)
%6,7,8)x,y,z (int,int,float) 
%9)diameter (int , each unit equals 1 micron) (pixels actually -A))
%10-)cell name 
%11 gfp sum
%
%For example:
%1, 1, 1, 1, -1, 305, 228, 10.2, 74, MS, 1602934, 0, 0, 0, , 
%
%USAGE:  M = readnuclei('t021-nuclei');
%        M = readnuclei('t021-nuclei', [1 3 4 5 6]); 
%OUTPUT : M (nX4) containing nuclei numeric index
%col1 - time point
%col2 - file specific index
%col3 - successor id(0 if the cell dies or divides, index of next timepoint)
%col456 - xyz
%7 -diameter
%8 gfp sum
%9 sucessor 1
%10 sucessor II
%11 predecessor
%READ FILE
fid = fopen(filename);

if fid < 0
     error(['could not open file: ' filename]);  
end

C = textscan(fid, '%f%d%d%f%f%f%f%f%f%s%f%f%f%f%s%f%f%f%f%f', 'delimiter', ',');
fclose(fid);

%PREPARE OUTPUT
if nargin == 1
     
%creating time vector from file name    
[pathstr, name] = fileparts(filename);
if (length(name) < 3)
    error('file name too short, could not extract the time point in which the file was taken')
end
t = str2double(name(2:4));
T = t*ones(length(C{1}),1);

%getting direct successors
temp=C{4};
temp((C{4} == -1) | (C{5} > 0)) = 0; %0 when the cell dies or divides
M = [T, C{1},temp, C{6},C{7},C{8},C{9},C{11},C{4},C{5},double(C{3}),C{12},C{13},C{19}];
cellNames = C{10};
forcedCellNames=C{15};
if  any( C{2} == -1)
    warning(['file : ' filename ' has invalid cells']);
end

valid_id =  find( C{2} == 1);  
M = M(valid_id,:);
cellNames = cellNames(valid_id);
else 
    error('col argument not implemented');
end

    
return;