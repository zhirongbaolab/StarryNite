# StarryNite
StarryNite is the current Matlab version of the StarryNite cell lineaging package, built initially for C. elegans but applied to a number of models.

To run StarryNite from source add the parent directory of the 3 source directories to the matlab path and run lineage_launcher_v2 ideally with distribution_lineaging as the working directory, this ensures that required .mat model files and jars for secondary tasks will be found.
If running from a compiled lineage_launcher_v2 the jars and .mat files in lineage launcher should be alongside the compiled executable. 

Code overview below:

Cell Detection
Directory:distribution_code.   Contains code for detection and segmentation of cells which can be used independently of lineaging.
  ProcessSequence process each image in a series loading the image and calling processVolume for each accumulating results in a data structure esequence that is ultimately output to a .mat file.  In ProcessVolume filtering is performed . slices are segmented and 3D maxima picked by createDiskSet. findOverlookedNuclei runs shape model to find unaccounted for nuclei. resolveConflicts examines all overlapping segmented nuclei to decide if they should be merged or remain split. 
The final internal data structure from this stage is described below and contains nuclear positions and a bunch of additional data about each nucleus detected much (but not necessarily all) of which is used during tracking. 
Note commandLineDriver is a historic top level interface for running segmentation from the command line, it should work, but is not currently used or supported.

Cell Lineaging
Directory:distribution_lineaging contains code which takes segmented nuclei and strings them together into a lineage
the main low level driver for tracking tracking_driver_new_classifier_based_version
which calls: 
initializeTrackingStructures to  initialize many of the fields used in tracking including calculating total gfp for nuclei and aspects of scoring the confidence that a cell is real
linkEasyCases to link easy cases and compute an explicit list of possible match candidates for each cell that is added to esequence and used subsequently (possible memory bottleneck?) 
greedyEndScore multiple times with different settings to create the tentative lineage (which systematically creates too many divisions)
greedydeleteFPbranches which resolves tentative bifurcations via classifier
Neither are used in the current UI called path but
trainingDriver is the top level driver for training bifurcation classifiers using existing detection result .mat files and  edited lineages (typically created by running the pipeline in naïve mode, or with an old model)
productionDriver is the top level driver for running normally, containing a subset of logic in the current top level driver for tracking and detection below. 

Driver UI
Directory: launcher_Interface  contains GUI used to interactively specify ROI for analysis and call cell detection and lineaging portions of pipeline.
Detect_track_driver_allmatlab in this directory is the most current top level driver for the current pipeline calling detection and tracking and generating all the files necessary to open the result in acetree.  This driver will work from the command line in theory (it just needs a few adjustments for the fact that the input arguments are all strings and not typed, these were made by bill but not applied here)

Data Structure
esequence (almost always referred to via this name as its passed around) is the main data structure used in both the detection and tracking code. It is a cell array of structures, with the length of the number of frames in the movie and a very large number of fields each of which is the length of the number of nuclei in that frame. Note that during detection the sizes of these fields changes. In this case some fields have the length of initial nuclei, and others have additional entries.  During tracking nuclei are not added or removed, false positives to be eliminated are indicated by an additional field (delete), and false negatives are handled via tracking between non consecutive timepoints, the actual interpolated points are generated while saving to acetree output.   
