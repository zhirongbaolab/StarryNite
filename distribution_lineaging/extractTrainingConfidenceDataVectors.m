function [ confidencedata] ...
    = extractTrainingConfidenceDataVectors( esequence,trackingparameters,embryonumbers_c,nucleidir,ROI,ROIxmin,ROIymin )
%compute all features for every link in tree and store for

    %match nuclei against reference having removed deleted nuclei
    %translate corrected successor and match back to full set
    esequence=matchLineageWithoutDeleted( esequence,trackingparameters,embryonumbers_c,nucleidir,ROI,ROIxmin,ROIymin );

    %rebuild answer key in order to get what matching thinks answer key
    %should be given deleted nuclei
    esequence=computeTrackingAnswerKey(esequence,trackingparameters.endtime);


nondiv=[];
div=[];
gap=[];
death=[];
nondivtrue=[];
divtrue=[];
gaptrue=[];
deathtrue=[];
confidencedata=[];
confidencedata.nondiv_ti=[];
confidencedata.div_ti=[];
confidencedata.gap_ti=[];
confidencedata.death_ti=[];
for t=1:trackingparameters.endtime-1
    stageinfo=mean(esequence{t}.selfdistance);
    if isinf(stageinfo)%single nucleus results in inf distance
        stageinfo=mean(esequence{t}.finaldiams);
    end
    for i=1:size(esequence{t}.finalpoints,1)
        if(~esequence{t}.delete(i)&max(esequence{t}.suc_time(i,:))<trackingparameters.endtime)
        if(esequence{t}.suc(i,1)~=-1)
            if (esequence{t}.suc(i,2)~=-1)
                %measure division
                divdata= computeDivisionConfidenceVector(esequence, t, i,trackingparameters);
                [match2,matchi,matchj ] = checkIfRealSuccessors( t,i,esequence{t}.suc_time(i,2),esequence{t}.suc(i,2),esequence );
                 [match1,matchi,matchj ] = checkIfRealSuccessors( t,i,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),esequence );
                 cdiv=match1&match2;
                 div=[div;divdata,stageinfo];
                 divtrue=[divtrue;cdiv];
                  confidencedata.div_ti=[confidencedata.div_ti;t,i];
            else
                if(esequence{t}.suc_time(i,1)~=t+1)
                    %measure simple nondiv
                    gapdata=     computeGapConfidenceVector(esequence, t, i,trackingparameters);
                    [match1,matchi,matchj ] = checkIfRealSuccessors( t,i,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),esequence );
                    gap=[gap;gapdata,stageinfo];
                    gaptrue=[gaptrue;match1&esequence{t}.correct_suc(i,2)==-1];
                     confidencedata.gap_ti=[confidencedata.gap_ti;t,i];
                else
                    %measure gap nondiv
                    nondivdata= computeNonDivisionConfidenceVector(esequence, t, i,trackingparameters);
                    [match1,matchi,matchj ] = checkIfRealSuccessors( t,i,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1),esequence );
                    nondiv=[nondiv;nondivdata,stageinfo];
                    nondivtrue=[nondivtrue;match1&esequence{t}.correct_suc(i,2)==-1];
                     confidencedata.nondiv_ti=[confidencedata.nondiv_ti;t,i];
                end
                
            end
        else
            %measure death
            deathdata=computeDeathConfidenceVector(esequence,t,i,trackingparameters);
            %death defined as predecessor, no sucessor and not in last 10
            %min of track, which is likely to be unrecovered FN
           % 'death better computed as matched nucleus has no sucessor'
            isdeath=esequence{t}.matches(i)~=-1&&...
                esequence{t}.corrected_sucessors(esequence{t}.matches(i),1)==-1;
           % isdeath=t<(trackingparameters.endtime-10)&~esequence{t}.FP(i)&&esequence{t}.correct_suc(i,1)==-1;
            death=[death;deathdata,stageinfo];
            deathtrue=[deathtrue;isdeath];
             confidencedata.death_ti=[confidencedata.death_ti;t,i];
        end
        end
    end
end
confidencedata.nondiv=nondiv;
confidencedata.div=div;
confidencedata.death=death;
confidencedata.gap=gap;
confidencedata.nondivtrue=nondivtrue;
confidencedata.divtrue=divtrue;
confidencedata.deathtrue=deathtrue;
confidencedata.gaptrue=gaptrue;

end
%{
%test code

%}
