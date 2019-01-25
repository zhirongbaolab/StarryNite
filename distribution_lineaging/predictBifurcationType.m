
function predicted_class = predictBifurcationType(...
    alldaughterdata,allforwarddata,allbackdata,d1length,d2length,...
    FNbackcand1lengths,FNbackcand2lengths,bestFNForwardLengthD1,...
    bestFNForwardLengthD2,bestFNBackCorrect,trackingparameters,bestIndex, count,forcemode)
%given all data about bifurcation broken into 3 blocks
%classify bifurcation with topology of bifurcation determining which blocks
%are present/relevant

%force classifier to not give up, picking the best option excluding the
%other class
if ~exist('forcemode')
    forcemode=false;
end

global computedclassificationvector;
global refclassificationvector;
global removed;
global simpleFNcorrect;

minsize=min(d1length,d2length);

smallcutoff=trackingparameters.smallcutoff;
%calculate which class this example falls in
FullyDivLooking=(d1length>=smallcutoff&d2length>=smallcutoff)&(~(FNbackcand1lengths>0)&~(FNbackcand2lengths>0));
FNDivLooking=(d1length>=smallcutoff&d2length>=smallcutoff)&(FNbackcand1lengths>0|FNbackcand2lengths>0);
small=minsize'<smallcutoff;
anysmalllacksforwardFN=(d1length)<smallcutoff&~(bestFNForwardLengthD1>0)|(d2length<smallcutoff&~(bestFNForwardLengthD2>0));
hasbackwardFN=(FNbackcand1lengths>0)|(FNbackcand2lengths>0);
FullyFPLooking=small&anysmalllacksforwardFN&~hasbackwardFN;
DirtyFPLooking=small&anysmalllacksforwardFN&hasbackwardFN;
DivFPLooking=  small&~anysmalllacksforwardFN&~hasbackwardFN;
trulyambigious=small&~anysmalllacksforwardFN&hasbackwardFN;

if(trulyambigious) %short, fn back, fn forward exist 
    data=[alldaughterdata(trackingparameters.bifurcationclassifier.daughterkeep),...
        allbackdata(trackingparameters.bifurcationclassifier.backkeep),...
        allforwarddata(trackingparameters.bifurcationclassifier.forwardkeep)];
 %   data(isinf(data))=10;
    predicted_class=predict(trackingparameters.bifurcationclassifier.ambigious,data,'HandleMissing','on');
 
    %force mode pick best non other class
    if(forcemode&predicted_class==0)
         posteriors_class=posterior(trackingparameters.bifurcationclassifier.ambigious,data,'HandleMissing','on');
            posteriors_class(1)=0;%null out other 
            [val,ind]=max(posteriors_class);
            classvals=[0,2,3];%hack i hate hard code this but seem to have a hard time extracting
            predicted_class=classvals(ind);
    end
    
 %no back, forward data is irrelevant because doesnt exist or long 
else if(FullyFPLooking|FullyDivLooking)
        data=alldaughterdata(trackingparameters.bifurcationclassifier.daughterkeep);
       % data(isinf(data))=10;
        predicted_class=predict(trackingparameters.bifurcationclassifier.fp_div,data,'HandleMissing','on');
        
        if(forcemode&predicted_class==0)
            posteriors_class=posterior(trackingparameters.bifurcationclassifier.fp_div,data,'HandleMissing','on');
            posteriors_class(1)=0;%null out other 
            [val,ind]=max(posteriors_class);
            classvals=[0,1,3];%hack i hate hard code this but seem to have a hard time extracting
            predicted_class=classvals(ind);
        end
        
        
        %back exists forward is irrelevant because doesnt exist or long
    else if(DirtyFPLooking|FNDivLooking)
            data=[alldaughterdata(trackingparameters.bifurcationclassifier.daughterkeep),...
                allbackdata(trackingparameters.bifurcationclassifier.backkeep)];
            %    data(isinf(data))=10;
            %           predicted_class=predict(trackingparameters.bifurcationclassifier.dirtyfp_fn,data);
            predicted_class=predict(trackingparameters.bifurcationclassifier.dirtyfp_fn,data,'HandleMissing','on');
            
            if(forcemode&predicted_class==0)
                posteriors_class=posterior(trackingparameters.bifurcationclassifier.dirtyfp_fn,data,'HandleMissing','on');
                posteriors_class(1)=0;%null out other 
            [val,ind]=max(posteriors_class);
            classvals=[0,1,2,3];%hack i hate hard code this but seem to have a hard time extracting
            predicted_class=classvals(ind);
            end
                  
        else
            if (DivFPLooking)%no back short and forward exists 
                data=[alldaughterdata(trackingparameters.bifurcationclassifier.daughterkeep),...
                    allforwarddata(trackingparameters.bifurcationclassifier.forwardkeep)];
         %       data(isinf(data))=10;
                predicted_class=predict(trackingparameters.bifurcationclassifier.divfp,data,'HandleMissing','on');
           
                if(forcemode&predicted_class==0)
                    posteriors_class=posterior(trackingparameters.bifurcationclassifier.divfp,data,'HandleMissing','on');
                   posteriors_class(1)=0;%null out other 
            [val,ind]=max(posteriors_class);
            classvals=[0,1,3];%hack i hate hard code this but seem to have a hard time extracting
            predicted_class=classvals(ind);

                end
                
            end
        end        
    end
end

%store computed class
computedclassificationvector=[computedclassificationvector,predicted_class];

if(forcemode&predicted_class==0)
         predicted_class=1;
end

%classifier sometimes fails in unexplored parts of space
if(isnan(predicted_class))
    'Nan returned as class label, unexpected or nan input input vector'
     if(forcemode)
         predicted_class=1;
     else
        predicted_class=0;
     end
end


%'short circut reall classiification with oracle'
%compute cheating answer online
%if( trackingparameters.trainingmode)
% as 2/19/2014 I added this to decouple storage of classification result from use of
% oracle in order to compute confusion matrix for test data
if(isfield(trackingparameters,'useclassifieroracle')&&trackingparameters.useclassifieroracle)
    %    
    cleanlybad=(removed(count-1,11)+removed(count-1,12)>=minsize');
    cleanlybad=cleanlybad&minsize<=smallcutoff;
   % slightlybad=(removed(count-1,11)+removed(count-1,12)>=minsize*.5');
   % slightlybad=slightlybad&minsize<=smallcutoff;
    realdiv=logical(removed(count-1,8));
    %    FNback=(simpleFNcorrect(count-1)==1)&~realdiv&~cleanlybad;
    FNback=(removed(count-1,23)|bestFNBackCorrect|(simpleFNcorrect(count-1)==1))&~realdiv&~cleanlybad;
    %fn back corect should be superset so was redundant above and can be taken out  
  %  FNback=(removed(count-1,23)|(simpleFNcorrect(count-1)==1))&~realdiv&~cleanlybad;
  
    % FNback=removed(count-1,23);%bestFNBackCorrect;%
    
   % other=(~cleanlybad&~FNback&~realdiv);
    
    
    predicted_class=0; %other=
    if(realdiv)
        predicted_class=1;   %div
    end
    %fn back
    if(~realdiv&FNback)
        predicted_class=2; %FN=
    end
    
    %fn back and correct
    %classtags(~realdiv&FNback&bestFNBackCorrect)=2; %FN=
    %fn back and backcorrect or other and fnbackcorrect
    %classtags(~realdiv&((FNback&bestFNBackCorrect)|(~realdiv&~FNback&~cleanlybad&bestFNBackCorrect)))=2; %FN=
    if (cleanlybad)
        predicted_class=3;
    end
    
end    %training mode end

%if predicted 2 and no fn option is other (only happens in training mode)
if(predicted_class==2&bestIndex==-1)
    predicted_class=0;
    'FN without match'
end

refclassificationvector=[refclassificationvector,predicted_class];


end

