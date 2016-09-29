%save file with paramters matching those in matlab param file, for
%starrynite
function outputSNParamFile(filename,xyres,zres,slices,starttime,endtime,nucsize)

fid=fopen(filename,'w');
   
  
fprintf(fid,'%s\r\n','# Sn parameter file automatically generated for loader version from matlab parameter file');
fprintf(fid,'%s\r\n','#');
fprintf(fid,'%s\r\n','# the start and end time points in the image series to be processed');
fprintf(fid,'%s %d\r\n','time_start',starttime);
fprintf(fid,'%s %d\r\n','time_end',endtime);


fprintf(fid,'%s\r\n','# the start and end planes in a stack to be processed');
fprintf(fid,'%s\r\n','plane_start 1');
fprintf(fid,'%s %d\r\n','plane_end',slices);

fprintf(fid,'%s\r\n','# image resolution (um)');
fprintf(fid,'%s %1.3f\r\n','xy_res',xyres);
fprintf(fid,'%s %1.2f\r\n','z_res',zres);

fprintf(fid,'%s\r\n','# time interval between time points (min)');
fprintf(fid,'%s\r\n','time_interval 1');

fprintf(fid,'%s\r\n','# expected diameter of nuclei for the first time point, measured in pixels');
fprintf(fid,'%s\r\n','# If your expected size is above 120, you should change the max_nuc_size to 200 (at the bottom of this file).');
fprintf(fid,'%s %d\r\n','nuc_size',nucsize); 

fprintf(fid,'%s\r\n','# expected diameter of the polar body, measured in pixels. Used to distinguish polar bodies from embryonic nuclei');
fprintf(fid,'%s %d\r\n','polar_size',nucsize/2);


fprintf(fid,'%s\r\n','# neighborhood size (for the low pass noise filter), use larger value for less noisy images');
fprintf(fid,'%s\r\n','neighborhood_size 15');


fprintf(fid,'%s\r\n','# various cutoffs for noise thresholding');

fprintf(fid,'%s\r\n','# the top fraction of the Gaussian noise distribution to be retained, i.e., 1 minus the noise_fraction is the fraction of the distribution to be masked.  The lower this value is, the more it masks');
fprintf(fid,'%s\r\n','noise_fraction 0.05');

fprintf(fid,'%s\r\n','# based on the GFP expression levels and image settings, we divide the series into three time zones and use different noise stringency.  The time_switches specify the boundaries of these time zones.  Measured in the number of cells present in the embryo.');
fprintf(fid,'%s\r\n','noise_time_switch1 55');
fprintf(fid,'%s\r\n','noise_time_switch2 180');

fprintf(fid,'%s\r\n','# additional noise cutoffs for finer control.  noise_factors are fudge factors applied to the noise thresholds.  nuc_density_cutoffs specify how much brighter a nucleus should be compared to the noise threshold.  Higher values of either may cause weak nuclei to be missed.  Three sets are provided for the three time zones.');
fprintf(fid,'%s\r\n','nuc_density_cutoff1 1.2');
fprintf(fid,'%s\r\n','noise_factor1 1.3');

fprintf(fid,'%s\r\n','nuc_density_cutoff2 1.6');
fprintf(fid,'%s\r\n','noise_factor2 1.3');

fprintf(fid,'%s\r\n','nuc_density_cutoff3 1.7');
fprintf(fid,'%s\r\n','noise_factor3 1.4');

fprintf(fid,'%s\r\n','# a cutoff to speed up the search for local maxima.  Higher value may cause weak nuclei to be missed.  Safest to set this to zero at the beginning of optimizing the noise cutoffs');
fprintf(fid,'%s\r\n','max_weight_cutoff 0.2');


fprintf(fid,'%s\r\n','# used to optimize nuclear size. (0,1].  Higher values make it more difficult to shrink or expand the spherical model of nuclei.  optimize this after optimizing the noise cutoffs. you might be able to further optimize the noise cutoff after optimizing these');
fprintf(fid,'%s\r\n','shrink_elastisity 0.6');
fprintf(fid,'%s\r\n','expand_elastisity 0.8');


fprintf(fid,'%s\r\n','# shortest cell cycle acceptable, in minutes');
fprintf(fid,'%s\r\n','minimal_cell_cycle 10');

fprintf(fid,'%s\r\n','# run time parameters that should not be changed lightly');

fprintf(fid,'%s\r\n','# used for relaxing the minimal movement algorithm.  Higher value means more nuclei at the previous time point would be considered as a potential match for the nucleus in question');
fprintf(fid,'%s\r\n','ambiguity_cutoff 1.2');


fprintf(fid,'%s\r\n','# maximal number of cells. Used to allocated memory to store the identified nuclei.  If the value is too low, the program will automatically adjust.  a value of 0 will cause the program to crash.');
fprintf(fid,'%s\r\n','cell_ct_limit 5000');


fprintf(fid,'%s\r\n','# used to calculate spehere models. If this value is smaller than the largest nucleus the program finds, it will cause trouble.  Measured in pixels');
fprintf(fid,'%s\r\n','max_nuc_size 150');


fprintf(fid,'%s\r\n','# window size and spatial range for the scanning box algorithm of nuclear identification, relative to nuc_size');
fprintf(fid,'%s\r\n','nuc_size_factor1 0.8');
fprintf(fid,'%s\r\n','nuc_size_factor2 0.8');
fprintf(fid,'%s\r\n','nuc_size_factor3 0.85');
fprintf(fid,'%s\r\n','nuc_size_factor4 0.75');


fprintf(fid,'%s\r\n','# run time flag.  changing graphic_output to 1 will lead the software to generate annotated images (with red circles to mark the spherical models of the identified nuclei).  ACeTree is a much better way to view the images and the annotation.  Mostly for debugging reasons.');
fprintf(fid,'%s\r\n','graphic_output 0');
fprintf(fid,'%s\n','');

 fclose(fid);