
function filtered=dogfilter(stack,sigma,znonuniformity)

%assure odd kernel so fft works right
%filter1o=fspecial3('gaussian',[floor(sigma/2)*2+1,floor(sigma/2)*2+1,max(3,floor(sigma/znonuniformity/2)*2+1)]);
%filter2o=fspecial3('gaussian',[floor(sigma*1.6/2)*2+1,floor(sigma*1.6/2)*2+1,max(3,floor(sigma*1.6/znonuniformity/2)*2+1)]);


%modified version decouples sigma and filter size
%so that at small sized filters doesnt blur too much in x bc 3 is smallest
%it can become
filter1o=fspecial3_mod('gaussian',[floor(sigma/2)*2+1,floor(sigma/2)*2+1,max(3,floor(sigma/znonuniformity/2)*2+1)],[sigma,sigma,sigma/znonuniformity]);
filter2o=fspecial3_mod('gaussian',[floor(sigma*1.6/2)*2+1,floor(sigma*1.6/2)*2+1,max(3,floor(sigma*1.6/znonuniformity/2)*2+1)],[sigma*1.6,sigma*1.6,sigma*1.6/znonuniformity]);

%Xfilts=imfilter(stack,filter1o,'replicate','same');
%Xfiltb=imfilter(stack,filter2o,'replicate','same');
%filtered=Xfilts-Xfiltb;

%Xfilts=imfilter(stack,fspecial3('gaussian',[sigma,sigma,sigma/znonuniformity]),'replicate','same');
%Xfiltb=imfilter(stack,fspecial3('gaussian',[sigma*1.6,sigma*1.6,sigma*1.6/znonuniformity]),'replicate','same');

%
Xfilts=conv3Dfreq(stack,filter1o);
filtered=conv3Dfreq(stack,filter2o);
filtered=Xfilts-filtered;

filtered=filtered./2; %fft is off by factor of 2
%[a,b]=memory;
%global parameters;
%parameters.memusage=max(parameters.memusage,a.MemUsedMATLAB);
