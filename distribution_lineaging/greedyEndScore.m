function esequence=greedyEndScore(esequence,trackingparameters)
%cereate links below a certain score threshold, used to create division and
%nondivision non easy 1-1 links, though it is capable of creating gap links
%if they are in the candidate set (currently they are not)

global answers;%store answer for post tunign
%greedyEndMatch(esequence,trackingparameters)
%  forward main cases for each end consider it dividing into up ahead
%  starts or continuing with or without a gap
nodes=0;

global alldatadiv;
global alldatanondiv;
global alldatadiv_taken;
global alldatanondiv_taken;

for t=trackingparameters.starttime:trackingparameters.endtime-1
   %hysteresis control loop, unused if no hysteresis;
    iterations=1;
  
    if(~trackingparameters.trackdiv&&isfield(trackingparameters,'hysteresis')&&trackingparameters.hysteresis)
        iterations=trackingparameters.hysteresisMaxsteps;
    end
    for iteration=1:iterations
          changed=false;
    for i=1:size(esequence{t}.finalpoints,1)
        %if this is an end or midpoint
        


        if(~isempty(esequence{t}.forwardcandidates{i})&&~esequence{t}.delete(i))
            nodes=nodes+1;
            
            % filter candidates to make sure still available
            candidates=esequence{t}.forwardcandidates{i}(:,1);
            candidates_t=esequence{t}.forwardcandidates{i}(:,2);
            
            %also filter to make sure has not been marked as deleted
            validcandidates=ones(size(candidates));
            for j=1:length(candidates)
                     
                %if hysteris is turned on checks 
                % that both  candidate and current cell are unlinked and below threshold
                %if so is a illegal match and is dropped
               hysteresisbad=isfield(trackingparameters,'hysteresis')&&...
                    trackingparameters.hysteresis&&...
                    esequence{t}.finalmaximas(i)<trackingparameters.hysteresis_intensityhigh&&...
                    esequence{t}.pred(i)==-1&&...
                    esequence{candidates_t(j)}.finalmaximas(candidates(j))<trackingparameters.hysteresis_intensityhigh&&...
                    esequence{candidates_t(j)}.suc(candidates(j),1)==-1;
                    
                
                if(hysteresisbad||...
                    esequence{candidates_t(j)}.pred(candidates(j))~=-1||...
                    esequence{candidates_t(j)}.delete(candidates(j)))
                    validcandidates(j)=0;
                end
            end
            candidates=candidates(logical(validcandidates));
            candidates_t=candidates_t(logical(validcandidates));
            
            %}
            
            scores_nondiv=[];
            
            
            %calculate lengths of all fragments involved in this
            endlength=-1;%traverse_back(esequence,t,i);
            %{
            startmerge=0;
            if(esequence{t}.fNN(i)>0)
                startmerge=length(esequence{t+1}.predecessor_suitors{esequence{t}.fNN(i)})>1 ...
                    |esequence{t+1}.bNN(esequence{t}.fNN(i))~=i;
            end
            
            endmerge=zeros(size(candidates,1),1);
            for k=1:size(candidates,1)
                if(esequence{candidates_t(k)}.bNN(candidates(k))~=-1)
                    endmerge(k)=length(esequence{candidates_t(k)-1}.sucessor_suitors{esequence{candidates_t(k)}.bNN(candidates(k))})>1|...
                        esequence{candidates_t(k)-1}.fNN(esequence{candidates_t(k)}.bNN(candidates(k)))~=candidates(k);
                    %lese remains
                end
            end
            %}
            %remove gaps that are too big, orr all of them
            good_nondivcandidates_ind=candidates_t-(t+1)<=trackingparameters.gapthresh;
            
            
            
            % good_nondivcandidates_ind=good_nondivcandidates_ind&gaptop;
            candidates=candidates(good_nondivcandidates_ind);
            candidates_t=candidates_t(good_nondivcandidates_ind);
            % end
            
            
            
            
            candidates_lengths=zeros(size(candidates));
            sucessor_length=-1;
            %{
            for k=1:length(candidates_lengths)
                candidates_lengths(k)=traverse_forward(esequence,candidates_t(k),candidates(k));
            end
            if(esequence{t}.suc(i,1)~=-1)
                sucessor_length=traverse_forward(esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1));
            end
            %}
            
            
            
            if(esequence{t}.suc(i,1)==-1&&trackingparameters.tracknondiv) %nondivision possible bc no links already
                %score cases via stored function handle
                
                   %{
                good_nondivcandidates_ind=candidates_t==t+1|(true);
                % good_nondivcandidates_ind=(candidates_lengths>2|endlength>1);
                
                candidates=candidates(good_nondivcandidates_ind);
                candidates_t=candidates_t(good_nondivcandidates_ind);
                candidates_lengths=candidates_lengths(good_nondivcandidates_ind);
                %}
                
                [scores_nondiv,splitscores,certainties_nondiv,nondivscorecomponents]=...
                    trackingparameters.nonDivCostFunction(esequence,i,t,candidates,candidates_t, trackingparameters);
                
                %  scores_nondiv=certainties_nondiv;
                
                if (trackingparameters.recordanswers)
                    %answer info
                    %anwer from cannonical answer
                    canswer_suc=esequence{t}.correct_suc(i,:);
                    canswer_suc_time=esequence{t}.correct_suc_time(i,:);
                    nondiv_choices_correct=zeros(size(scores_nondiv));
                    top2=zeros(size(scores_nondiv));
                    %
                    for j=1:length(nondiv_choices_correct)
                        %tests if it is correct in sense of being nondiv and
                        %correct
                        %{
                    nondiv_choices_correct(j)=canswer_suc(2)==-1&&...
                        canswer_suc(1)==candidates(j)&&...
                        canswer_suc_time(1)==candidates_t(j);
                        %}
                        
                        %tests in sense of being part of real lineage
                        %(returns true for half of a division)
                        nondiv_choices_correct(j)=...
                            (canswer_suc(1)==candidates(j)&&...
                            canswer_suc_time(1)==candidates_t(j))|...
                            (canswer_suc(2)==candidates(j)&&...
                            canswer_suc_time(2)==candidates_t(j));
                        top2(j)=length(esequence{candidates_t(j)}.predecessor_suitors{candidates(j)});
                        
                    end
                    %end ID, correct,  scores( s1, c1 c2),size1 size2, gap size
                    %top1, top2
                    top1=length(esequence{t}.sucessor_suitors{i})*ones(size(splitscores,1),1);
                    
                    if(~isempty(splitscores))
                        alldatanondiv=[alldatanondiv;ones(size(splitscores,1),1)*nodes,nondiv_choices_correct,splitscores,endlength*ones(size(splitscores,1),1),candidates_lengths,candidates_t-t,top1,top2,nondivscorecomponents];
                    end
                else
                    if(~isempty(splitscores))
                        alldatanondiv=[alldatanondiv;ones(size(splitscores,1),1)*nodes,ones(size(splitscores,1),1)*-1,splitscores];
                    end
                    
                end
                
                
                
            end
            divisionlengths=[];
            %score divisions
            endtop=length(esequence{t}.sucessor_suitors{i});
            divcandidates_t=[];
            divcandidates=[];
            if(trackingparameters.trackdiv)
                if(esequence{t}.suc(i,1)==-1) %divisions are pairs possible bc no links already
                    if(trackingparameters.do_enddiv)
                        %expand candidates to all pairs
                        
                        
                        for j=1:length(candidates)-1
                            for k=j+1:length(candidates)
                                divcandidates_t=[divcandidates_t;candidates_t(j),candidates_t(k)];
                                divcandidates=[divcandidates;candidates(j),candidates(k)];
                                divisionlengths=[divisionlengths;endlength,candidates_lengths(j),candidates_lengths(k)];
                            end
                        end
                        divtypes=ones(size(divcandidates,1),1);
                    end
                else %divisions are existing link and all candidates
                    
                    divcandidates=[candidates,esequence{t}.suc(i,1)*ones(size(candidates))];
                    divcandidates_t=[candidates_t,esequence{t}.suc_time(i,1)*ones(size(candidates_t))];
                    divtypes=2*ones(size(divcandidates,1),1);
                    divisionlengths=[divisionlengths;endlength*ones(size(candidates_lengths)),candidates_lengths,sucessor_length*ones(size(candidates_lengths))];
                    
                end
            end
            
            
          
            %score divisions
            
            [scores_div,splitscores_div,certainties_div]=trackingparameters.DivCostFunction(esequence,i,t,divcandidates,divcandidates_t, trackingparameters);
            %scores_div=certainties_div;
            
            %data on all candidate scorings
            
            if(~isempty(divcandidates)&trackingparameters.recordanswers)
                
                dtop=[];
                for cand=1:size(divcandidates,1)
                    dtop=[dtop;endtop,length(esequence{divcandidates_t(cand,1)}.predecessor_suitors{divcandidates(cand,1)}),length(esequence{divcandidates_t(cand,2)}.predecessor_suitors{divcandidates(cand,2)})];
                end
                
                canswer_suc=esequence{t}.correct_suc(i,:);
                canswer_suc_time=esequence{t}.correct_suc_time(i,:);
                div_choices_correct=zeros(size(scores_div));
                for j=1:length(div_choices_correct)
                    onecorrect=(canswer_suc(1)==divcandidates(j,1)&&canswer_suc_time(1)==divcandidates_t(j,1))...
                        |(canswer_suc(1)==divcandidates(j,2)&&canswer_suc_time(1)==divcandidates_t(j,2));
                    twocorrect=(canswer_suc(2)==divcandidates(j,1)&&canswer_suc_time(2)==divcandidates_t(j,1))...
                        |(canswer_suc(2)==divcandidates(j,2)&&canswer_suc_time(2)==divcandidates_t(j,2));
                    
                    div_choices_correct(j)=onecorrect&twocorrect;
                end
                
                
                
                %end ID, sucessor exists, correct, scores( agreement s1 s2, c1 c2 c3),size1 size2 size3 gap size
                alldatadiv=[alldatadiv;ones(size(divtypes))*nodes,divtypes,div_choices_correct,splitscores_div,divisionlengths,divcandidates_t-t,dtop];
                
            end
            
            %take min
            [minscore_div,divimin]=min(scores_div);
            [minscore_nondiv,imin]=min(scores_nondiv);
            
            if (trackingparameters.recordanswers)
                %answer info
                %anwer from cannonical answer
                canswer_suc=esequence{t}.correct_suc(i,:);
                canswer_suc_time=esequence{t}.correct_suc_time(i,:);
                
                %is answer in candidate set at all
                answerpresent=0;
  
                
                minscore_nondivs=minscore_nondiv;
                minscore_divs=minscore_div;
                if (isempty(imin))
                    minscore_nondivs=-1;
                end
                if (isempty(divimin))
                    minscore_divs=-1;
                end
                %is min score in each category the answer
                if(canswer_suc(2)==-1)
                    divanswerright=0;
                    if (~isempty(imin))
                        nondivanswerright=min(canswer_suc==[candidates(imin,:),-1]&canswer_suc_time==[candidates_t(imin,:),-1]);
                    else
                        nondivanswerright=-1;
                    end
                else% division in reality
                    nondivanswerright=0;
                    if (~isempty(imin))
                        nondivanswerright=(canswer_suc(1)==candidates(imin)&canswer_suc_time(1)==candidates_t(imin))...
                            |(canswer_suc(2)==candidates(imin)&canswer_suc_time(2)==candidates_t(imin));
                        
                    end
                    if (~isempty(divimin))
                        divanswerright=min((canswer_suc==divcandidates(divimin,:)))&min(canswer_suc_time==divcandidates_t(divimin,:));
                        flip=divcandidates(divimin,:);
                        flip(1)=divcandidates(divimin,2); flip(2)=divcandidates(divimin,1);
                        flipt(1)=divcandidates_t(divimin,2); flipt(2)=divcandidates_t(divimin,1);
                        divanswerright=divanswerright|(min((canswer_suc==flip))&min(canswer_suc_time==flipt));
                        
                    else
                        
                        divanswerright=-1;
                    end
                    
                end
                
                
            end %end if should record answers
            
            %actually do tracking
            %overwrite scores used to pick min with confidences in their
            %correctness for thresholding/ div/nondiv selection
            %{
            if(trackingparameters.lookup)
                if(~isempty(minscore_div))
                    minscore_div=certainties_div(divimin);
                    %do lookup into calibrated max certainty table
                    minscore_div=histogram_lookup(trackingparameters.model.final_div_bins,trackingparameters.model.final_div_fraction_wrong,minscore_div);
                    minscore_divs=minscore_div;
                end
                if(~isempty(minscore_nondiv))
                    minscore_nondiv=certainties_nondiv(imin);
                    minscore_nondiv=histogram_lookup(trackingparameters.model.final_nondiv_bins,trackingparameters.model.final_nondiv_fraction_wrong,minscore_nondiv);
                    minscore_nondivs=minscore_nondiv;%make stored one lookedup
                end
            end
            %}
            minscore=[];
            div=true;
            nondiv=true;
            %would apply thresholds here if wanted
            if(isempty(minscore_div)||minscore_div>trackingparameters.endscorethresh_div)
                div=false;
            end
            if(isempty(minscore_nondiv)||minscore_nondiv>trackingparameters.endscorethresh_nondiv)
                nondiv=false;
            end
            if(nondiv&&div)
                if(minscore_div<minscore_nondiv)
                    % if(trackingparameters.divalpha_constant+minscore_div*trackingparameters.divalpha<minscore_nondiv)
                    div=true;
                    nondiv=false;
                else
                    div=false;
                    nondiv=true;
                end
            end
            
            
            
            if(div||nondiv)
                %record answer only when you deal with node
                %min score div, min score nondiv,
                %divright,nondivright,answerpresent (always 0),answer, answertime, i,t
                if(trackingparameters.recordanswers)
                    %score of correct answer
                    if(canswer_suc(2)==-1)
                        if(canswer_suc(1)~=-1)
                            [correctscore,splitscoresc]=trackingparameters.nonDivCostFunction(esequence,i,t,canswer_suc(1),canswer_suc_time(1), trackingparameters);
                        else
                            correctscore=-1;
                        end
                    else
                        correctscore=trackingparameters.DivCostFunction(esequence,i,t,canswer_suc,canswer_suc_time, trackingparameters);
                    end
                    
                    %anwswer stuff
                    if(div)
                        canswer=[divcandidates(divimin,1:2),divcandidates_t(divimin,1:2)];
                        sucFN=min(esequence{divcandidates_t(divimin,1)}.correct_suc(divcandidates(divimin,1),1),...
                            esequence{divcandidates_t(divimin,2)}.correct_suc(divcandidates(divimin,2),1));
                        
                    else
                        canswer=[candidates(imin),-1,candidates_t(imin),-1];
                        sucFN=esequence{candidates_t(imin)}.correct_suc(candidates(imin),1);
                        
                    end
                    answers=[answers;...
                        minscore_divs,minscore_nondivs,divanswerright,nondivanswerright,answerpresent,canswer_suc,canswer_suc_time,i,t,sucFN,correctscore,canswer];
                end
                
                if(trackingparameters.dotrack) %actually do tracking instead of gathering case info
                    %actual linking back on
                    if(div)%div best score
                        changed=true;
                         
                        esequence{t}.suc(i,1:2)=divcandidates(divimin,1:2);
                        esequence{t}.suc_time(i,1:2)=divcandidates_t(divimin,1:2);
                        esequence{divcandidates_t(divimin,1)}.pred(divcandidates(divimin,1))=i;
                        esequence{divcandidates_t(divimin,2)}.pred(divcandidates(divimin,2))=i;
                        esequence{divcandidates_t(divimin,1)}.pred_time(divcandidates(divimin,1))=t;
                        esequence{divcandidates_t(divimin,2)}.pred_time(divcandidates(divimin,2))=t;
                        
                        %remove linked candidates from foward candidates lists of
                        %each of its back candidates should go here rather than
                        %redundantly filtering list
                        predcandidates=esequence{divcandidates_t(divimin,1)}.backcandidates{divcandidates(divimin,1)};
                        
                        for j=1:size(predcandidates,1) %remove linked sucessor from each of its previous
                            
                            forwardcandidates=esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)};
                            if (~isempty(forwardcandidates)) %could be empty bc answered and ccleared toally
                                equalentry=forwardcandidates(:,1)==divcandidates(divimin,1)&...
                                    forwardcandidates(:,2)==divcandidates_t(divimin,1);
                                esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)}=...
                                    forwardcandidates(~equalentry,:);
                            end
                        end
                        
                        predcandidates=esequence{divcandidates_t(divimin,2)}.backcandidates{divcandidates(divimin,2)};
                        for j=1:size(predcandidates,1) %remove linked sucessor from each of its previous
                            
                            forwardcandidates=esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)};
                            if (~isempty(forwardcandidates)) %could be empty bc answered and ccleared toally
                                equalentry=forwardcandidates(:,1)==divcandidates(divimin,2)&...
                                    forwardcandidates(:,2)==divcandidates_t(divimin,2);
                                esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)}=...
                                    forwardcandidates(~equalentry,:);
                            end
                        end
                        
                        %if we make a division mark node as completed
                        esequence{t}.forwardcandidates{i}=[];
                        %}
                    else %nondiv best score
                     changed=true;
                        
                        %  'creating nondivision'
                        esequence{t}.suc(i,1)=candidates(imin);
                        esequence{t}.suc_time(i,1)=candidates_t(imin);
                        esequence{candidates_t(imin)}.pred(candidates(imin))=i;
                        esequence{candidates_t(imin)}.pred_time(candidates(imin))=t;
                        %remove linked candidates from foward candidates lists of
                        %each of its back candidates should go here rather than
                        %redundantly filtering list
                        
                        predcandidates=esequence{candidates_t(imin)}.backcandidates{candidates(imin)};
                        for j=1:size(predcandidates,1) %remove linked sucessor from each of its previous
                            forwardcandidates=esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)};
                            if (~isempty(forwardcandidates)) %could be empty bc answered and ccleared toally
                                equalentry=forwardcandidates(:,1)==candidates(imin)&...
                                    forwardcandidates(:,2)==candidates_t(imin);
                                esequence{predcandidates(j,2)}.forwardcandidates{predcandidates(j,1)}=...
                                    forwardcandidates(~equalentry,:);
                            end
                        end
                        %}
                    end
                end
                
                
            end %best score is better than threshold
            
        end %does it have candidates
    end
    if (~changed)
        break
    end
    end%iterations
end


end

