%linearly interpolates from highthresh in upper right corner of volume
%to lowthresh in lower right corner 
function twodmax=DLSMAdaptiveMaximaFilter(volsize,twodmax,maximavals,numcells)
%{
highthresh=getParameter('a_highthresh',numcells);
lowthresh=getParameter('a_lowthresh',numcells);
[x,y,z]=ind2sub(volsize,twodmax);
threshvals=(y/volsize(2)).^2;%fraction of way to far right
threshvals=threshvals.*(z/volsize(3)).^2;
threshold=threshvals*highthresh+(1-threshvals)*lowthresh;
twodmax=twodmax(maximavals>threshold);
%}

highthreshx=getParameter('a_highthreshx',numcells);
lowthreshx=getParameter('a_lowthreshx',numcells);
highthreshz=getParameter('a_highthreshz',numcells);
lowthreshz=getParameter('a_lowthreshz',numcells);
[x,y,z]=ind2sub(volsize,twodmax);
threshvalsx=(y/volsize(2)).^2;%fraction of way to far right
threshvalsz=(z/volsize(3)).^2;

threshold=.5*(threshvalsx*highthreshx+(1-threshvalsx)*lowthreshx);
threshold=threshold+(.5*(threshvalsz*highthreshz+(1-threshvalsz)*lowthreshz));
twodmax=twodmax(maximavals>threshold);
