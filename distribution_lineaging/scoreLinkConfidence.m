function [ esequence] ...
    = scoreLinkConfidence(esequence,trackingparameters)
%compute all features for every link in tree computeconfidence and store

%initialize 
for t=1:trackingparameters.endtime-1
esequence{t}.linkconfidences=zeros(size(esequence{t}.finalpoints,1),3);
%esequence{t}.linkconfidences=zeros(size(esequence{t}.finalpoints,1),1);
end

for t=1:trackingparameters.endtime-1
        stageinfo=mean(esequence{t}.selfdistance);

    for i=1:size(esequence{t}.finalpoints,1)
        if(~esequence{t}.delete(i)&max(esequence{t}.suc_time(i,:))<=trackingparameters.endtime)
        if(esequence{t}.suc(i,1)~=-1)
            if (esequence{t}.suc(i,2)~=-1)
                %measure division
                divdata= [computeDivisionConfidenceVector(esequence, t, i,trackingparameters),stageinfo];
               confidence=posterior(trackingparameters.linkconfidencemodel.div,...
                    divdata(:,trackingparameters.linkconfidencemodel.divkeep));
                 confidence(isnan(confidence))=1;%nan is presumed wrong
                 confidence(isinf(confidence))=1;%nan is presumed wrong
                %lookup actual empirical probability of incorrect draw from
                % bin
               confidence=trackingparameters.linkconfidencemodel.divlookup...
                   (round(confidence(1)*(length(trackingparameters.linkconfidencemodel.bins)-1))+1) ;
          % confidence=confidence(1);
            esequence{t}.linkconfidences(i,1)=1-confidence;
            
            %bif side 1
                                nondivdata= [computeNonDivisionConfidenceVector(esequence, t, i,1,trackingparameters),stageinfo];
                    confidence=posterior(trackingparameters.linkconfidencemodel.nondiv,...
                        nondivdata(:,trackingparameters.linkconfidencemodel.nondivkeep));
                    confidence(isnan(confidence))=1;%nan is presumed wrong
                    confidence(isinf(confidence))=1;%nan is presumed wrong
                    %for some reason I made this histogram the opposite way
                    %around
                    if((confidence(2))>trackingparameters.linkconfidencemodel.nondivlowerlimit&...
                            (confidence(2))<trackingparameters.linkconfidencemodel.nondivupperlimit)
                        confidence=trackingparameters.linkconfidencemodel.nondivlookup...
                            (round((confidence(2))*(length(trackingparameters.linkconfidencemodel.nondivlookup)-1))+1);
                    else
                        if((confidence(2))<=trackingparameters.linkconfidencemodel.nondivlowerlimit)
                            confidence=trackingparameters.linkconfidencemodel.nondivlowerlimitp;
                        else
                            confidence=trackingparameters.linkconfidencemodel.nondivupperlimitp;
                        end
                    end
                      esequence{t}.linkconfidences(i,2)=confidence;
                      %bif side 2
                                                      nondivdata= [computeNonDivisionConfidenceVector(esequence, t, i,2,trackingparameters),stageinfo];
                    confidence=posterior(trackingparameters.linkconfidencemodel.nondiv,...
                        nondivdata(:,trackingparameters.linkconfidencemodel.nondivkeep));
                    confidence(isnan(confidence))=1;%nan is presumed wrong
                    confidence(isinf(confidence))=1;%nan is presumed wrong
                    %for some reason I made this histogram the opposite way
                    %around
                    if((confidence(2))>trackingparameters.linkconfidencemodel.nondivlowerlimit&...
                            (confidence(2))<trackingparameters.linkconfidencemodel.nondivupperlimit)
                        confidence=trackingparameters.linkconfidencemodel.nondivlookup...
                            (round((confidence(2))*(length(trackingparameters.linkconfidencemodel.nondivlookup)-1))+1);
                    else
                        if((confidence(2))<=trackingparameters.linkconfidencemodel.nondivlowerlimit)
                            confidence=trackingparameters.linkconfidencemodel.nondivlowerlimitp;
                        else
                            confidence=trackingparameters.linkconfidencemodel.nondivupperlimitp;
                        end
                    end
                      esequence{t}.linkconfidences(i,3)=confidence;
                    
       
            else
                if(esequence{t}.suc_time(i,1)~=t+1)
                    %measure gap  nondiv
                    gapdata=[computeGapConfidenceVector(esequence, t, i,trackingparameters),stageinfo];
                    confidence=posterior(trackingparameters.linkconfidencemodel.gap,...
                        gapdata(:,trackingparameters.linkconfidencemodel.gapkeep));
                     confidence(isnan(confidence))=1;%nan is presumed wrong
                        confidence(isinf(confidence))=1;%nan is presumed wrong
            %          confidence=confidence(1);
                    confidence=trackingparameters.linkconfidencemodel.gaplookup...
                        (round(confidence(1)*(length(trackingparameters.linkconfidencemodel.bins)-1))+1) ;
                  esequence{t}.linkconfidences(i,1)=1-confidence;
         
                else
                    %measure nondiv
                    nondivdata= [computeNonDivisionConfidenceVector(esequence, t, i,1,trackingparameters),stageinfo];
                    confidence=posterior(trackingparameters.linkconfidencemodel.nondiv,...
                        nondivdata(:,trackingparameters.linkconfidencemodel.nondivkeep));
                    confidence(isnan(confidence))=1;%nan is presumed wrong
                    confidence(isinf(confidence))=1;%nan is presumed wrong
                    %for some reason I made this histogram the opposite way
                    %around
                    if((confidence(2))>trackingparameters.linkconfidencemodel.nondivlowerlimit&...
                            (confidence(2))<trackingparameters.linkconfidencemodel.nondivupperlimit)
                        confidence=trackingparameters.linkconfidencemodel.nondivlookup...
                            (round((confidence(2))*(length(trackingparameters.linkconfidencemodel.nondivlookup)-1))+1);
                    else
                        if((confidence(2))<=trackingparameters.linkconfidencemodel.nondivlowerlimit)
                            confidence=trackingparameters.linkconfidencemodel.nondivlowerlimitp;
                        else
                            confidence=trackingparameters.linkconfidencemodel.nondivupperlimitp;
                        end
                    end
                      esequence{t}.linkconfidences(i,1)=confidence;
                    
                    %    confidence=trackingparameters.linkconfidencemodel.nondivlookup...
                    %         (round(confidence(1)*(length(trackingparameters.linkconfidencemodel.bins)-1))+1) ;
                    %      confidence=confidence(1);
                    %      esequence{t}.linkconfidences(i)=1-confidence;
                end
                
            end
        else
            %measure death
            deathdata=[computeDeathConfidenceVector(esequence,t,i,trackingparameters),stageinfo];
              confidence=posterior(trackingparameters.linkconfidencemodel.death,...
                       deathdata(:,trackingparameters.linkconfidencemodel.deathkeep));
                    confidence(isnan(confidence))=1;%nan is presumed wrong
                   confidence=trackingparameters.linkconfidencemodel.deathlookup...
                       (round(confidence(1)*(length(trackingparameters.linkconfidencemodel.bins)-1))+1) ;
                  esequence{t}.linkconfidences(i,1)=confidence;
                   %confidence=confidence(1);
                %    esequence{t}.linkconfidences(i)=1-confidence;
        end
        end
    end
end


end
%{
%test code

%}
