%confidence model training below here
%assemble data for link confidence
alllinkconfidencedata=allbifurcationinfo(1).linkconfidencedata;
for lin=2:length(allbifurcationinfo)
    alllinkconfidencedata.nondiv=[alllinkconfidencedata.nondiv;allbifurcationinfo(lin).linkconfidencedata.nondiv];
    alllinkconfidencedata.div=[alllinkconfidencedata.div;allbifurcationinfo(lin).linkconfidencedata.div];
    alllinkconfidencedata.death=[alllinkconfidencedata.death;allbifurcationinfo(lin).linkconfidencedata.death];
    alllinkconfidencedata.gap=[alllinkconfidencedata.gap;allbifurcationinfo(lin).linkconfidencedata.gap];
     alllinkconfidencedata.nondivtrue=[alllinkconfidencedata.nondivtrue;allbifurcationinfo(lin).linkconfidencedata.nondivtrue];
    alllinkconfidencedata.divtrue=[alllinkconfidencedata.divtrue;allbifurcationinfo(lin).linkconfidencedata.divtrue];
    alllinkconfidencedata.deathtrue=[alllinkconfidencedata.deathtrue;allbifurcationinfo(lin).linkconfidencedata.deathtrue];
    alllinkconfidencedata.gaptrue=[alllinkconfidencedata.gaptrue;allbifurcationinfo(lin).linkconfidencedata.gaptrue];
end
model=trainConfidenceClassifier(alllinkconfidencedata);
trackingparameters.linkconfidencemodel=model;
bins=linspace(0,1,10);
trackingparameters.linkconfidencemodel.bins=bins;
%test model


%div
confidence_div=posterior(trackingparameters.linkconfidencemodel.div,...
    alllinkconfidencedata.div(:,trackingparameters.linkconfidencemodel.divkeep));
distance_div=(alllinkconfidencedata.div(:,5).^2+alllinkconfidencedata.div(:,6).^2).^.5;
cases=alllinkconfidencedata.divtrue;

scores=confidence_div;
%scores=distance_div;


acases=alllinkconfidencedata.divtrue;
ascores=confidence_div;
ascoresdis=distance_div;

ghc=histc((scores(logical(cases),1)),bins);
bhc=histc((scores(~logical(cases),1)),bins);
ratio=bhc./(ghc+bhc);

fill=ratio(find(~isnan(ratio),1));
for i=1:length(bins)
    if (isnan(ratio(i)))
     ratio(i)=fill;
    else
        fill=ratio(i);
    end
end

trackingparameters.linkconfidencemodel.divlookup=ratio;


%nondiv
confidence_nondiv=posterior(trackingparameters.linkconfidencemodel.nondiv,...
    alllinkconfidencedata.nondiv(:,trackingparameters.linkconfidencemodel.nondivkeep));
distance_nondiv=(alllinkconfidencedata.nondiv(:,3).^2+alllinkconfidencedata.nondiv(:,4).^2).^.5;
cases=alllinkconfidencedata.nondivtrue;

scores=confidence_nondiv;
%scores=distance_nondiv;

acases=[acases;alllinkconfidencedata.nondivtrue];

ascores=[ascores;confidence_nondiv];
ascoresdis=[ascoresdis;distance_nondiv];

%bins=linspace(0,.00000000000000000000001,10);
bins=linspace(0,1,10);

ghc=histc((scores(logical(cases),2)),bins);
bhc=histc((scores(~logical(cases),2)),bins);
ratio=bhc./(ghc+bhc)


ghc=histc((scores(logical(cases),1)),bins);
bhc=histc((scores(~logical(cases),1)),bins);
ratio=bhc./(ghc+bhc);

fill=ratio(find(~isnan(ratio),1));
for i=1:length(bins)
    if (isnan(ratio(i)))
     ratio(i)=fill;
    else
        fill=ratio(i);
    end
end

trackingparameters.linkconfidencemodel.nondivlookup=ratio;

%gap
confidence_gap=posterior(trackingparameters.linkconfidencemodel.gap,...
    alllinkconfidencedata.gap(:,trackingparameters.linkconfidencemodel.gapkeep));
distance_gap=(alllinkconfidencedata.gap(:,3).^2+alllinkconfidencedata.gap(:,4).^2).^.5;
cases=alllinkconfidencedata.gaptrue;

scores=confidence_gap;
%scores=distance_gap;


acases=[acases;alllinkconfidencedata.gaptrue];
ascores=[ascores;confidence_gap];
ascoresdis=[ascoresdis;distance_gap];


ghc=histc((scores(logical(cases),1)),bins);
bhc=histc((scores(~logical(cases),1)),bins);
ratio=bhc./(ghc+bhc);

fill=ratio(find(~isnan(ratio),1));
for i=1:length(bins)
    if (isnan(ratio(i)))
     ratio(i)=fill;
    else
        fill=ratio(i);
    end
end

trackingparameters.linkconfidencemodel.gaplookup=ratio;

%death
confidence_death=posterior(trackingparameters.linkconfidencemodel.death,...
    alllinkconfidencedata.death(:,trackingparameters.linkconfidencemodel.deathkeep));
cases=alllinkconfidencedata.deathtrue;

scores=confidence_death;
%scores=distance_gap;
ghc=histc((scores(logical(cases),1)),bins);
bhc=histc((scores(~logical(cases),1)),bins);
ratio=bhc./(ghc+bhc);

fill=ratio(find(~isnan(ratio),1));
for i=1:length(bins)
    if (isnan(ratio(i)))
     ratio(i)=fill;
    else
        fill=ratio(i);
    end
end

trackingparameters.linkconfidencemodel.deathlookup=ratio;
return

%ROC all
figure
[sortascores,IX]=sort(ascores(:,1));
sortacases=acases(IX);
badc=[];
goodc=[];
for i=1:length(sortacases);
    goodc=[goodc;sum(sortacases(1:i))];
    badc=[badc;sum(~sortacases(1:i))];
end    
plot( badc./sum(~sortacases),goodc./sum(sortacases));
hold on

[sortascores,IX]=sort(ascoresdis(:,1));
sortacases=acases(IX);
badc=[];
goodc=[];
for i=1:5:length(sortacases);
    goodc=[goodc;sum(sortacases(1:i))];
    badc=[badc;sum(~sortacases(1:i))];
end    
plot( badc./sum(~sortacases),goodc./sum(sortacases),'g');




bins=linspace(0,1,10);
ghc=histc((ascores(logical(acases),1)),bins);
bhc=histc((ascores(~logical(acases),1)),bins);
figure 
plot(bins,bhc./(ghc+bhc),'r');
hold on
plot(bins,ghc./(ghc+bhc),'b');

figure
plot(ghc,'b');
figure
plot(bhc,'r');
%roc

bins=linspace(0,1,100);
ghc=histc((scores(logical(cases),1)),bins);
bhc=histc((scores(~logical(cases),1)),bins);
goodc=zeros(size(bins));
badc=zeros(size(bins));
for i=1:length(bins)
goodc(i)=sum(ghc(1:i));
badc(i)=sum(bhc(1:i));
end
plot( badc./sum(bhc),goodc./sum(ghc));

%score vs fraction correct  
figure 
plot(bins,bhc./ghc);

%combine over all classes