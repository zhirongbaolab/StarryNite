function [ cells,embinfouned] = loadcells_unnamed(directory,endtime,anisotropy,loadcon )
%loads cells in directory 
%returning names expressions and a dummy list of matched names
%(for historical reasons of running in debugging with answer key mode)
%this version is not dependent on name hashing for collecting timepoints
%but uses the links

[embinfouned,errorsuned ]...
    = loadEmbryo_unzipped( directory,endtime);

for time=1:endtime
    %if (time<length(embinfouned))&& ~isempty(embinfouned(time))
    embinfouned(time).expression=embinfouned(time).celldata(:,12)-embinfouned(time).celldata(:,14);
    embinfouned(time).pred=-1*ones(size(embinfouned(time).expression));
    
    if(loadcon)
        if (exist([directory,'confidence\t',num2str(time,'%03d'),'-con.mat']))
            confidence=load([directory,'confidence\t',num2str(time,'%03d'),'-con.mat']);
            confidence_t=ones(size(embinfouned(time).expression));
            if(~isempty(confidence.con))
                confidence_t(1:size(confidence.con,1),:)=confidence.con(:,1);
            end
            confidence=confidence_t;
        else
            confidence=ones(size(embinfouned(time).expression));
        end
    else
        confidence=ones(size(embinfouned(time).expression));
    end
    embinfouned(time).confidence=confidence;
    %end
end

for time=1:endtime
    % if (time<length(embinfouned))&& ~isempty(embinfouned(time))
      for i=1:size(embinfouned(time).finalpoints,1)
            if(embinfouned(time).suc(i,1)~=-1)
                 embinfouned(time+1).pred(embinfouned(time).suc(i,1))=i;
            end
            if(embinfouned(time).suc(i,2)~=-1)
                 embinfouned(time+1).pred(embinfouned(time).suc(i,2))=i;
            end
      end
     % end
end

cells={};
cellcount=1;
for t=1:endtime
   %  if (t<length(embinfouned))&& ~isempty(embinfouned(t))
    for i=1:size(embinfouned(t).finalpoints,1)
        %parse new cell if from div or start
        if(cellcount==366)
        'odd'
        end
        if(embinfouned(t).pred(i)==-1||embinfouned(t-1).suc(embinfouned(t).pred(i),2)~=-1)
            continueparse=true;
            count=[];
            count.exp=[];
            count.con=[];
            count.pos=[];
            count.matchnames={};
            count.divides=false;
            count.birthposition=embinfouned(t).finalpoints(i,:);
            count.name=embinfouned(t).names{i};
            currenti=i;
            currentt=t;
            while(continueparse)
                if( currentt~=endtime&& embinfouned(currentt).suc(currenti,2)~=-1&& embinfouned(currentt).suc(currenti,1)~=-1)
                    count.divides=true;
                    count.sucnames={embinfouned(currentt+1).names{ embinfouned(currentt).suc(currenti,1)},...
                        embinfouned(currentt+1).names{embinfouned(currentt).suc(currenti,2)}};
                    continueparse=false;
                    count.endtime=currentt;
                else if (currentt==endtime||embinfouned(currentt).suc(currenti,1)==-1)
                        continueparse=false;
                        count.endtime=currentt;
                    end
                end
                count.pos=[count.pos;embinfouned(currentt).finalpoints(currenti,:)];
                count.exp=[count.exp;embinfouned(currentt).expression(currenti)];
                count.con=[count.con;embinfouned(currentt).confidence(currenti)];
                count.matchnames={count.matchnames{:},' '};
                currenti=embinfouned(currentt).suc(currenti,1);
                currentt=currentt+1;
            end
            cells{cellcount}=count;
            cellcount=cellcount+1;
        end
   % end
      end
end