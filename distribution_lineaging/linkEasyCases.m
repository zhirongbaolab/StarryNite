function esequence=linkEasyCases(esequence,trackingparameters)
%link easy 1-1 cases (currently mutual non conflicting)
global answers_1;
global answers;
answers_1=[];
answers=[];

% solves  frame to frame easy cases
%initialize a data structure for things that think nucleus should be its
%predecessor all passed back in esequence

for t=trackingparameters.starttime:trackingparameters.endtime
esequence{t}.sucessor_suitors=cell(1,size(esequence{t}.finalpoints,1));
esequence{t}.sucessor_suitors_score=cell(1,size(esequence{t}.finalpoints,1));
esequence{t}.predecessor_suitors=cell(1,size(esequence{t}.finalpoints,1));
esequence{t}.predecessor_suitors_score=cell(1,size(esequence{t}.finalpoints,1));
esequence{t}.sucessor_suitors_splitscore=cell(1,size(esequence{t}.finalpoints,1));
esequence{t}.fNN=-1*ones(1,size(esequence{t}.finalpoints,1));
esequence{t}.bNN=-1*ones(1,size(esequence{t}.finalpoints,1));
end

%build a list of all claimants to a nuclesu as predecessor based on nn of
%score
for t=trackingparameters.endtime:-1:trackingparameters.starttime+1
    if(~(isempty(esequence{t-1}.finalpoints)|isempty(esequence{t}.finalpoints)))
    distances=distance_anisotropic(esequence{t-1}.finalpoints',esequence{t}.finalpoints',trackingparameters.anisotropyvector);

    for i=1:size(esequence{t}.finalpoints,1)
        
        candidates=linspace(1,size(esequence{t-1}.finalpoints,1),size(esequence{t-1}.finalpoints,1));
    %    candidates_t=t-1*ones(size(candidates));
        
        
        %score cases via stored function handle
     %   scores=[];
     %   splitscores=[];
       % for s=1:length(candidates)
      %      [score_c,splitscores_c]=trackingparameters.nonDivCostFunction(esequence,candidates(s),t-1,i,t, trackingparameters);
      %      scores=[scores,score_c];
      %      splitscores=[splitscores;splitscores_c];
      %  end
      scores=distances(:,i);
      splitscores=scores;
      
        [minscore,imin]=min(scores);
        %concatenate best candidate into long sucessor vector
        if(~isempty(imin))
            esequence{t}.bNN(i)=candidates(imin);
            esequence{t-1}.sucessor_suitors{candidates(imin)}=[esequence{t-1}.sucessor_suitors{candidates(imin)};i];
            esequence{t-1}.sucessor_suitors_score{candidates(imin)}=[esequence{t-1}.sucessor_suitors_score{candidates(imin)};minscore];
            esequence{t-1}.sucessor_suitors_splitscore{candidates(imin)}=[esequence{t-1}.sucessor_suitors_splitscore{candidates(imin)};splitscores(imin,:)];
            
        end
        
    end
    end
end

%build same list in opposite direction
for t=trackingparameters.starttime:trackingparameters.endtime-1
   if(~(isempty(esequence{t}.finalpoints)|isempty(esequence{t+1}.finalpoints)))
  
    distances=distance_anisotropic(esequence{t}.finalpoints',esequence{t+1}.finalpoints',trackingparameters.anisotropyvector);

    for i=1:size(esequence{t}.finalpoints,1)
       % candidates=esequence{t}.forwardcandidates{i};
         candidates=linspace(1,size(esequence{t+1}.finalpoints,1),size(esequence{t+1}.finalpoints,1));
      candidates_t=t+1*ones(size(candidates));
         %score cases via stored function handle
        %scores=[];
        %splitscores=[];
        %for s=1:length(candidates)
         %   [scores,splitscores]=trackingparameters.nonDivCostFunction(esequence,i,t,candidates,candidates_t, trackingparameters);
        %    scores=[scores,score_c];
        %     splitscores=[splitscores;splitscores_c];
        %end
         scores=distances(i,:);
        [minscore,imin]=min(scores);

        %concatenate best candidate into long sucessor vector
        if(~isempty(imin))
            esequence{t}.fNN(i)=candidates(imin);
        esequence{t+1}.predecessor_suitors{candidates(imin)}=[esequence{t+1}.predecessor_suitors{candidates(imin)};i];
        esequence{t+1}.predecessor_suitors_score{candidates(imin)}=[esequence{t+1}.predecessor_suitors_score{candidates(imin)};minscore];
        end
    end
   end
end


for t=trackingparameters.starttime:trackingparameters.endtime-1
    for i=1:size(esequence{t}.finalpoints,1)
        target=[];
        
        
        mutual=false;
        for k=1:length(esequence{t}.sucessor_suitors{i})
            for l=1:length(esequence{t+1}.predecessor_suitors{esequence{t}.sucessor_suitors{i}(k)})
            if(esequence{t+1}.predecessor_suitors{esequence{t}.sucessor_suitors{i}(k)}(l)==i)
                mutual=true;
                tcand=k;
                fcand=l;
            end
            end
        end
        %mutual and no forward or back conflict
        
        if(trackingparameters.conflictfilter)
            goodlink=(mutual&&...
                length(esequence{t+1}.predecessor_suitors{esequence{t}.sucessor_suitors{i}(tcand)})==1&&...
                length(esequence{t}.sucessor_suitors{i})==1);
        else
            goodlink=mutual;
        end
        
        if(goodlink)
            
        %mutual no forward conflict and not just after a forward conflict
        % if(mutual&&length(esequence{t+1}.predecessor_suitors{esequence{t}.sucessor_suitors{i}(tcand)})==1 ...
         %      &&  length(esequence{t}.predecessor_suitors{i})<=1 )
    
           %        &&length(esequence{t}.sucessor_suitors{i})<=2)
%mutual
          % if(mutual)
      
            target=esequence{t}.sucessor_suitors{i}(tcand);
            targetscore=esequence{t}.sucessor_suitors_score{i}(tcand);
            
            splitscore=esequence{t}.sucessor_suitors_splitscore{i}(tcand,:);
            
            safefactor=trackingparameters.safefactor;
            %safe distanceverride only makes sense when using distance as
            %cost function**
            mdistance=distanceCostFunction(esequence,i,t,target,t+1, trackingparameters)*mean(esequence{t}.selfdistance);
            safedistance=esequence{t}.selfdistance(i)/safefactor;%dist to nn/2
            safedistance=min(safedistance,esequence{t+1}.selfdistance(target)/safefactor);
            
            if(trackingparameters.safefilter)
                if(~isempty(imin))
                    
                    if(mdistance>safedistance) %wipe out all candidates if best is unsafe
                        target=[];
                        targetscore=[];
                    end
                end
                
            end
             hysteresisbad=true;
            if(~isempty(target))
            hysteresisbad=isfield(trackingparameters,'hysteresis')&&...
                trackingparameters.hysteresis&&...
                esequence{t}.finalmaximas(i)<trackingparameters.hysteresis_intensityhigh&&...
                esequence{t}.pred(i)==-1&&...
                esequence{t+1}.finalmaximas(target)<trackingparameters.hysteresis_intensityhigh&&...
                esequence{t+1}.suc(target,1)==-1;
            end
            
            
            %if(~isempty(target)&&targetscore<trackingparameters.nondivscorethreshold(t))%&...
            if(~hysteresisbad&&~isempty(target)&&targetscore<trackingparameters.forwardcutoff(t))%&...
                
                esequence{t}.suc(i,1)=target;
                esequence{t}.suc_time(i,1)=t+1;
                esequence{t+1}.pred(target)=i;
                esequence{t+1}.pred_time(target)=t;
                %remove successor from possible sucessor list of everything
                %on its possible predecessor list;
%{
                predcandidates=esequence{t+1}.backcandidates{target};
                for j=1:size(predcandidates,1) %remove linked sucessor from each of its previous
                    forwardcandidates=esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)};
                    equalentry=forwardcandidates(:,1)==target&...
                        forwardcandidates(:,2)==t+1;
                    esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)}=...
                        forwardcandidates(~equalentry,:);
                    
                end
                esequence{t+1}.backcandidates{target}=[];
  %}              
                
                if( trackingparameters.recordanswers)
                    %answer calculation
                    canswer_suc=esequence{t}.correct_suc(i,:);%the true successors to chosen predecessor
                    canswer_suc_time=esequence{t}.correct_suc_time(i,:);%time
                    answerpresent=0;
                    correct=0;
                    if(~isempty(target))
                        correct=((canswer_suc(1)==target&&canswer_suc_time(1)==t+1)||(canswer_suc(2)==target&&canswer_suc_time(2)==t+1));
                    end
                    
                    % next time point, best score, best score is
                    %correct, local thresold, is the predecessor linked to anything in
                    %answer?,is current linked to anything forward in answer, time i
                    %candidate(i) secondscore, distance safedistance,
                    %splitscore, fp start, fp end, fn end
                    if(isempty(targetscore))
                        'why would this happen?'
                        answers_1=[answers_1;-1,-1,correct,trackingparameters.nondivscorethreshold(t),-1,-1,esequence{t}.correct_suc(i,:),...
                            t,i,-1,-1,-1,-1,-1,-1,-1,-1,-1];
                        %    answers_1=[answers_1;-1,-1,correct,trackingparameters.nondivscorethreshold(t),-1,-1,esequence{t}.correct_suc(i,:),...
                        %      t,i,-1,-1,-1,-1,-1];
                    else
                        % secondscore=min(scores(scores>minscore));
                        % if(isempty(secondscore))
                        %     secondscore=-1;
                        % end
                        answers_1=[answers_1;canswer_suc_time(1)==t+1|canswer_suc_time(2)==t+1,targetscore,correct,trackingparameters.nondivscorethreshold(t),...
                            canswer_suc,esequence{t}.correct_suc(i,:),t,i,target,-1,mdistance,safedistance,splitscore,esequence{t}.FP(i),esequence{t+1}.FP(target),esequence{t+1}.correct_pred_time(target)>0&&esequence{t+1}.correct_pred_time(target)~=t];
                    end
                end
                
            end
        end
  
        
        
    end
end
end

