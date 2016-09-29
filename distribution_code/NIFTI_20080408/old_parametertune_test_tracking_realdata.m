clear all
%load data files
% convert them into tracking data

allresults={};
c1=0;

for var1=2.5:.5:5
    c1=c1+1;
    c2=0;
    for var2=2.5:.5:5
        c2=c2+1;
        
testfile='L:\duz\project\Imaging\BV24_NewJournalV_20091008\ZD_BV24_JournalV_1_s1_tracktest_sn21_fullmatlabresult.mat';
load (testfile);
time=300;
startt=0
movieInfo=[];
%struct(time,1);
%pos=rand(10,3)*500;
%amp=rand(10,1)*2;
clear backup_esequence;
global parameters;

for i=startt+1:time
   [integratedGFP,area]=integrateGFP(esequence{i},parameters);

   %amp=esequence{i}.finalmaximas;  
   amp=single(integratedGFP);
%   pos=esequence{i}.finalaveragepoints(1:length(amp),:);
   pos=esequence{i}.finalpoints(1:length(amp),:);
   pos(:,3)=pos(:,3)*anisotropy;
%{
movieInfo(i-startt).xCoord=[pos(:,1),zeros(size(pos(:,1)))];
movieInfo(i-startt).yCoord=[pos(:,2),zeros(size(pos(:,2)))];
movieInfo(i-startt).zCoord=[pos(:,3),zeros(size(pos(:,3)))];
movieInfo(i-startt).amp=[amp(:,1),zeros(size(amp(:,1)))];
  %}
   if (~ROI)
       ROIxmin=0;
       ROIymin=0;
   end
        movieInfo(i-startt).xCoord=[pos(:,1)+ROIxmin,zeros(size(pos(:,1)))];
        movieInfo(i-startt).yCoord=[pos(:,2)+ROIymin,zeros(size(pos(:,2)))];
        movieInfo(i-startt).zCoord=[pos(:,3),zeros(size(pos(:,3)))];
        movieInfo(i-startt).amp=[amp,zeros(size(amp))];
end

test={};
for i=1:length(esequence)
    test{i}.finaldiams=esequence{i}.finaldiams;
end
clear amp
clear pos
clear integratedGFP
clear area
clear esequence
esequence=test;




%% Cost functions

%compiler note to include these files that are not explicitly called
%#function costMatLinearMotionLink costMatLinearMotionCloseGaps kalmanResMemLM kalmanInitLinearMotion kalmanGainLinearMotion

%Frame-to-frame linking
costMatrices(1).funcName = 'costMatLinearMotionLink';

%Gap closing, merging and splitting
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';

%--------------------------------------------------------------------------

%% Kalman filter functions

%Memory reservation
kalmanFunctions.reserveMem = 'kalmanResMemLM';

%Filter initialization
kalmanFunctions.initialize = 'kalmanInitLinearMotion';

%Gain calculation based on linking history
kalmanFunctions.calcGain = 'kalmanGainLinearMotion';

%--------------------------------------------------------------------------

%% General tracking parameters

%Gap closing time window
gapCloseParam.timeWindow = 5;%10;

%Flag for merging and splitting
gapCloseParam.mergeSplit = 1;

%Minimum track segment length used in the gap closing, merging and
%splitting step
gapCloseParam.minTrackLen = 2;

%--------------------------------------------------------------------------

%% Cost function specific parameters: Frame-to-frame linking

%Flag for linear motion
parameters.linearMotion = 0;

%Search radius lower limit
parameters.minSearchRadius = 2;%2

%Search radius upper limit
parameters.maxSearchRadius = 50;%100;%5

%Standard deviation multiplication factor
parameters.brownStdMult =var1;% 2.4;%5;%8

%Flag for using local density in search radius estimation
parameters.useLocalDensity = 1;

%Number of past frames used in nearest neighbor calculation
parameters.nnWindow = gapCloseParam.timeWindow;

%Store parameters for function call
costMatrices(1).parameters = parameters;
clear parameters

%--------------------------------------------------------------------------

%% Cost cunction specific parameters: Gap closing, merging and splitting

%Same parameters as for the frame-to-frame linking cost function
parameters.linearMotion = costMatrices(1).parameters.linearMotion;
parameters.useLocalDensity = costMatrices(1).parameters.useLocalDensity;
parameters.maxSearchRadius = costMatrices(1).parameters.maxSearchRadius;
parameters.minSearchRadius = costMatrices(1).parameters.minSearchRadius;

%*** I changed thiis brownian std /2 before 1 is default .75 for test 1 .5
%conservative %2
parameters.brownStdMult = var2*costMatrices(1).parameters.brownStdMult*ones(gapCloseParam.timeWindow,1);
parameters.nnWindow = costMatrices(1).parameters.nnWindow;

%Gap length (frames) at which f(gap) (in search radius definition) reaches its
%plateau
parameters.timeReachConfB = 2;

%Amplitude ratio lower and upper limits
parameters.ampRatioLimit = [0.1 8];%[.25, 4] [0.5 4];

%Minimum length (frames) for track segment analysis
parameters.lenForClassify = 4;%5

%Standard deviation multiplication factor along preferred direction of
%motion
parameters.linStdMult = var2*ones(gapCloseParam.timeWindow,1);
%parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1);
%Gap length (frames) at which f'(gap) (in definition of search radius
%parallel to preferred direction of motion) reaches its plateau
parameters.timeReachConfL = gapCloseParam.timeWindow;

%Maximum angle between the directions of motion of two linear track
%segments that are allowed to get linked
parameters.maxAngleVV = 90;%45;

%Store parameters for function call
costMatrices(2).parameters = parameters;
clear parameters

%--------------------------------------------------------------------------




%track

scriptTrackGeneral_noparam;
%plotTracks2D(tracksFinal,[], '1', [],  1, 1)

%evaluate
mkdir('nuclei');

saveNucleiFiles(tracksFinal,esequence,time,'nuclei',anisotropy)

scripted_evaluateStarrynite_with_tracking
rmdir('nuclei','s');


allresults{c1,c2}=results;
    end
end
