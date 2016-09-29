function FFT1=conv3Dfreq(cprocl,vker)

%pad volume cprocl with zeros ov vker/2 width all around so that weird
%wrapping artifact does not influence output
%{
kernelwidth=floor(size(vker)/2);
cprocl_e=zeros(size(cprocl)+kernelwidth*2);
centerker=floor(size(cprocl)/2);
centerpix=floor(size(cprocl_e)/2);
embed1im=centerpix(1)-centerker(1)+[1:size(cprocl,1)];
embed2im=centerpix(2)-centerker(2)+[1:size(cprocl,2)];
embed3im=centerpix(3)-centerker(3)+[1:size(cprocl,3)];
cprocl_e(embed1im,embed2im,embed3im)=cprocl;
cprocl=cprocl_e;
%}
s=size(vker);
cprocl=padarray(cprocl,[floor(s(1)/2),floor(s(2)/2),floor(s(3)/2)],'replicate');
%cprocl=padarray(cprocl,[0,0,floor(s(3)/2)],'replicate');

%cprocl=double(cprocl);
smoothcell=zeros(size(cprocl));
centerpix=floor(size(cprocl)/2);
centerker=floor(size(vker)/2);

embed1i=centerpix(1)-centerker(1)+[1:size(vker,1)];
embed2i=centerpix(2)-centerker(2)+[1:size(vker,2)];
embed3i=centerpix(3)-centerker(3)+[1:size(vker,3)];

smoothcell(embed1i,embed2i,embed3i)=vker;

% disp('matched filter correlating observed volume');
% Correlation through FFT needs:
% a. Finding FFT(S) and FFT(C)

%pad with zero to kernel size so never get wrapping?
%FFT1=fftn(cprocl,size(cprocl)+2*floor(size(vker)/2));
%FFT2=fftn(smoothcell,size(cprocl)+2*floor(size(vker)/2));

FFT1=fftn(cprocl);
FFT2=fftn(smoothcell);

% b. ComplexArray=FFT(S) * Conj{FFT(C)}
%CPXARR2=FFT1.*conj(FFT2);
% c. IFFT{ComplexArray}
%MFout=2*abs(fftshift(ifftn(CPXARR2)));

% c. IFFT{ComplexArray}
%MFout=2*(fftshift(ifftn(FFT1.*conj(FFT2))));
FFT1=2*abs(fftshift(ifftn(FFT1.*conj(FFT2))));
s=floor(size(vker)/2);
s2=size(FFT1);

FFT1=FFT1(s(1)+1:s2(1)-s(1),s(2)+1:s2(2)-s(2),s(3)+1:s2(3)-s(3));
%FFT1=FFT1(:,:,s(3)+1:s2(3)-s(3));

%MFout=Mfout(embed1im,embed2im,embed3im); %crop away padding
%{

tic
Xfilts=imfilter(Xorig,fspecial3('gaussian',[sigma,sigma,sigma/anisotropy]),'replicate','same');
Xfiltb=imfilter(Xorig,fspecial3('gaussian',[sigma*1.6,sigma*1.6,sigma*1.6/anisotropy]),'replicate','same');

filtered=Xfilts-Xfiltb;
time=toc
tic
Xfilts2=conv3Dfreq(Xorig,fspecial3('gaussian',[sigma,sigma,sigma/anisotropy]));
Xfiltb2=conv3Dfreq(Xorig,fspecial3('gaussian',[sigma*1.6,sigma*1.6,sigma*1.6/anisotropy]));
filtered2=Xfilts2-Xfiltb2;
time=toc

%}