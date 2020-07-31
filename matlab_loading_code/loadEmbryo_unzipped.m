function [embinfo,errors ]...
    = loadEmbryo_unzipped(file,endtime)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


errors=[];


% examine all divisions and see if A cell is more anterior based on long
% axis of variance and axis label in auxinfo file
%templocation='temp_unzip\';
%{
%unzip zipfile to temp file
if ~exist(zipfile,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(zipfile,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end
%}
embinfo=[];
allpoints=[];

%load entire embryo
templocation=file;
for t=1:endtime
    
    %load t current
    nuclei=[templocation,'nuclei\t',num2str(t,'%03d'),'-nuclei'];
    if ~exist(nuclei)
           nuclei=[templocation,'nuclei\t',num2str(t,'%03d'),'-nuclei.txt'];
          if ~exist(nuclei)  
            errors.nucleiFiles=['expected nuclei file missing time: ',nuclei];
            return
          end
    end
    [celldata,cellnames,forcednames]=readnuclei_allfields(nuclei);
    embinfo(t).celldata=celldata;
    embinfo(t).cellnames=cellnames;
    p1_sucessors=celldata(:,9:10);
    celllocations=celldata(:,4:6);%pull nuclei from labeled data
    ABaindex=find(strcmp('ABa',cellnames));
    P2index=find(strcmp('P2',cellnames));
    if(~isempty(ABaindex)&&~isempty(P2index))
        ABapos=celllocations(ABaindex,:);
        P2pos=celllocations(P2index,:);
    end
    
    
    
    
    if(t<endtime)
        %load t+1
        nuclei=[templocation,'nuclei\t',num2str(t+1,'%03d'),'-nuclei'];
        if ~exist(nuclei)
              nuclei=[templocation,'nuclei\t',num2str(t+1,'%03d'),'-nuclei.txt'];
               if ~exist(nuclei)
                errors.nucleiFiles=['expected nuclei file missing time: ',nuclei];
                return
               end
        end
        [celldata_c2,cellnames2]=readnuclei(nuclei);
        indiciesp2=celldata_c2(:,2);
        
        %translate acetree sucessor indicies to array indicies
        s=size(p1_sucessors);
        p1_sucessors_t=[];
        for j=1:s(1);
            suc1=-1;
            suc2=-1;
            if(p1_sucessors(j,1)~=-1)
                suc1=find(indiciesp2==p1_sucessors(j,1));
                if isempty(suc1)
                    suc1=-1;
                    'invalid nucleus pointed to by valid nucleus, should not happen in new files'
                end
            end
            if(p1_sucessors(j,2)~=-1)
                suc2=find(indiciesp2==p1_sucessors(j,2));
                if isempty(suc2)
                    suc2=-1;
                    'invalid nucleus pointed to by valid nucleus, should not happen in new files'
                end
            end
            p1_sucessors_t=[p1_sucessors_t;suc1,suc2];
        end
    else
        p1_sucessors_t=p1_sucessors;% dont try and translate if endtime
    end
    embinfo(t).finalpoints=celllocations;
    embinfo(t).suc=p1_sucessors_t;
    embinfo(t).names=cellnames;
    embinfo(t).forcednames=forcednames;
      

end

end

