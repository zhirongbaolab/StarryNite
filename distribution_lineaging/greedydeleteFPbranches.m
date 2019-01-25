function esequence=greedydeleteFPbranches(esequence,trackingparameters)
%walk through tentative naive lineage,
% optionally eliminate truly isolated points that cant form bifurcations
% judge each bifuration and based on judgment of its origin (unknown, fn,
% fp, division) modify linkage of tree around bifurcation

global removed;%store answer for post tunign
global FNtype;
global simpleFNcorrect;
global removedi;
global BifurcationMeasures;
global candidateInfo;
global bestCandidateInfo;
global endIndicies;
global bestEndCandidateInfo;
global confidenceData;
global splitFNMatchScore;
global classround;
count=1;
%first pass removing truly isolated points which are left dangling because
%they have no candidates
if(trackingparameters.deleteisolated)
    for t=trackingparameters.starttime:trackingparameters.endtime-1
        
        for i=1:size(esequence{t}.finalpoints,1)
            if(~esequence{t}.delete(i))
                
                %if this is a candidateless fragment
                if(esequence{t}.pred(i)==-1)
                    if(esequence{t}.suc(i,1)~=-1)
                        d1length=traverse_forward(esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1));
                        d1cur=esequence{t}.suc(i,1);
                        tcur1=esequence{t}.suc_time(i,1);
                        d1sumconfidence=0;
                        d1FP=0;
                        for j=1:d1length
                            % d1sumconfidence=d1sumconfidence+esequence{tcur1}.confidences(d1cur);
                            if(isfield(esequence{tcur1},'FP'))
                                d1FP=d1FP+esequence{tcur1}.FP(d1cur);
                            end
                            d1curbk=d1cur;
                            d1cur=esequence{tcur1}.suc(d1cur,1);
                            tcur1=esequence{tcur1}.suc_time(d1curbk,1);
                        end
                        numcells=length(esequence{t}.finaldiams);
                        removedi=[removedi;d1length,d1sumconfidence,d1FP,numcells];
                        
                        delete=((d1length<trackingparameters.FPsizethresh&numcells<=trackingparameters.earlythresh)...
                            | d1length<=trackingparameters.FPsizethreshsmall);
                        %(log(d1sumconfidence/d1length)>trackingparameters.deletethresh));
                        %delete=false;
                        if(delete)
                            
                            d1cur=esequence{t}.suc(i,1);
                            tcur1=esequence{t}.suc_time(i,1);
                            for j=1:d1length %mark short branch for deletion
                                esequence{tcur1}.delete(d1cur)=1;
                                d1curbk=d1cur;
                                d1cur=esequence{tcur1}.suc(d1cur,1);
                                tcur1=esequence{tcur1}.suc_time(d1curbk,1);
                            end
                            if (esequence{t}.suc(i,2)==-1)
                                esequence{t}.suc_time(i,1)=-1;
                                esequence{t}.suc(i,1)=-1;
                                esequence{t}.delete(i)=1;
                            else %if says to delete branch 1 as orphan and is a branch 2 swap them  and dont delete curent
                                esequence{t}.suc_time(i,1)=esequence{t}.suc_time(i,2);
                                esequence{t}.suc(i,1)=esequence{t}.suc(i,2);
                                % esequence{esequence{t}.suc_time(i,2)}.pred(esequence{t}.suc(i,2))=i;
                                % esequence{esequence{t}.suc_time(i,2)}.pred(esequence{t}.suc(i,2))=t;
                                esequence{t}.suc_time(i,2)=-1;
                                esequence{t}.suc(i,2)=-1;
                            end
                            
                        end
                        
                    else%if it is totally isolated delete it
                        esequence{t}.delete(i)=1;
                        numcells=length(esequence{t}.finaldiams);
                        if (trackingparameters.recordanswers)
                            removedi=[removedi;1,esequence{t}.confidences(i),esequence{t}.FP(i),numcells];
                        end
                    end
                end
            end
        end
    end
end%end check whether to delete isolated

%iterate over whole lineage looking for bifurcations when found
%gather data
%classify
%take appropriate action
dFNf=0;
dFNb=0;
for t=trackingparameters.starttime:trackingparameters.endtime-1
    
    %control logic for whether to force bifurcations to do soemthing used
    %primarily for forcing an answer during early timepoints
    if(isfield(trackingparameters,'useforce')&&trackingparameters.useforce&&t<trackingparameters.forceendtime)
        forcestatus=true;
    else
        forcestatus=false;
    end
    
    for i=1:size(esequence{t}.finalpoints,1)
        %if this is a division
        if(esequence{t}.suc(i,2)~=-1&~esequence{t}.delete(i))
            
            %assemble data about bifurcation
            [d1cand,d2cand,d1candt,d2candt,d1length,d2length,...
                bestFNBackCorrect,bestmatchings,FNbackcand1lengths,FNbackcand2lengths,...
                bestFNForwardLengthD1,bestFNForwardLengthD2,....
                bestmatchingsplayerstart,bestmatchingplayerend,bestdaughter,bestIndex,...
                splitscores_div,alldaughterdata,allforwarddata,allbackdata,count] ....
                = assembleBifurcationData( esequence,t,i,trackingparameters,count );
           
       %     if(t==79&i==17)
       %     'blah'
       %     end
              %pick whether to call the single, or multiple statistical model
            %classifier based on whether it is a single or multi model
            %passed in by parameter file
            if (isfield(trackingparameters.bifurcationclassifier,'ambigious'))
        
            predicted_class = predictBifurcationType(...
                alldaughterdata,allforwarddata,allbackdata,d1length,d2length,...
                FNbackcand1lengths,FNbackcand2lengths,bestFNForwardLengthD1,...
                bestFNForwardLengthD2,bestFNBackCorrect,trackingparameters,bestIndex, count,forcestatus);
            else
                 predicted_class = predictBifurcationTypeSinglemodel(...
                alldaughterdata,allforwarddata,allbackdata,d1length,d2length,...
                FNbackcand1lengths,FNbackcand2lengths,bestFNForwardLengthD1,...
                bestFNForwardLengthD2,bestFNBackCorrect,trackingparameters,bestIndex, count,forcestatus);
            end
            minsize=min(d1length,d2length);
            classround=[classround;1];
            if (trackingparameters.dotrack)
                %process bifurcation based on predicted class
                if(predicted_class==0)
                    %other
                    FNtype=[FNtype;-1, -1];
                    [esequence,count] = processOtherBifurcation( esequence,t,i,splitscores_div,trackingparameters,count );
                else
                    if(predicted_class==3)
                        %FP
                        FNtype=[FNtype;-1, -1];
                        esequence = processFPBifurcation( esequence,t,i,minsize,d1length );
                    else
                        if (predicted_class==2)
                            %FN
                            
                            [esequence,numplayers,endplayerssameasdiv ]= processFNBifurcation...
                                ( esequence,t,i,bestmatchings, bestmatchingsplayerstart,...
                                bestmatchingplayerend, bestdaughter,bestIndex,d1cand,d2cand,...
                                d1candt,d2candt,trackingparameters);
                            
                            %  matchinfo
                            FNtype=[FNtype;numplayers, endplayerssameasdiv];
                        else
                            %division
                            FNtype=[FNtype;-1, -1];
                        end %dont need to do anthing for class=1 division
                    end
                end
            end
            
            
        end
    end
end

end

