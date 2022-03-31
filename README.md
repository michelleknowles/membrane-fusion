# membrane-fusion
Protocol: Find and Measure Single Fusion Events

The purpose of this code is to automate the identification of fusion events in fluorescene microscopy data (time series) that give rise to an increase in intensity, typically due to the unquenching of a pH sensitive fluorophore. We have used it with pHluorin and pHuji tagged to CD63, a multi-vesicular endosome marker. It is not a GUI and requires some knowledge of Matlab. 
Approach: The general workflow is described in a figure named "workflow.pdf". We analyze one movie at a time, but that could be streamlined. Our program uses other freely available code. 

Required:

•	Software: Matlab and Matlab Toolbox “Image Processing”

•	Files (here): stkwrite, stkread, CellMaskMaker3, FusionDataCompile, FusionEventFinderV5, MiniStackFusionDataProcessV9, maskMaker3, LoadMovie, FusionMiniMaker8, CellIntensityInTime

•	Other Files (not ours): pkfind, bpass, shared here: https://site.physics.georgetown.edu/matlab/

•	Add all of the above to your path. 


1.	MAKE A STORAGE FOLDER: We use a four-digit identifier unique for each movie (3242, in this example). Go to this folder.
2.	ASSIGN FILENAME: Assign the movie that shows fusion and the move that shows a protein of interest or other color channel. This is my fusion filename ‘C1-3242-5 A549 CD63pHluorin Syx4myc594 37C 50ms 300g 0ND transformed croppedR.stk’. You need the stk or tif extension here. Make this a shorter name to handle in calls: 
  >> fusion_name = 'C1-3242-5 A549 CD63pHluorin Syx4myc594 37C 50ms 300g 0ND transformed croppedR.stk';
  >> protein_name = 'C2-3242-5 A549 CD63pHluorin Syx4myc594 37C 50ms 300g 0ND transformed croppedL.stk';
3.	LOAD MOVIES:
  >> fusion_movie = LoadMovie(fusion_name);
  
  >> protein_movie = LoadMovie(protein_name);
4.	MAKE A CELL MASK AND MEASURE BACKGROUND, AVERAGE CELL INTENSITIES: Make a mask of the cell to obtain the out of cell background and the average cell intensity. Manually adjust the threshold so that it looks like the cell is isolated from the background. The threshold depends on cell intensity. Not essential for finding fusion events but used to subtract a background. This puts two variables into the workspace, writes out two files to the active working directory and plots the images below. This makes a mask from the first image. Note this works for our data because the plasma membrane is fluorescent with our fusion marker. You can use a different channel if you have a marker that is solely localized to a vesicle (NPY/tPA). 
  >> [cellMask, backgroundMask, avgcellintensity] = CellMaskMaker3(fusion_movie, 0.3, 'cellmask_name.stk');

If you need the cell intensity for every frame in time (for when photobleaching is an issue). Note that this uses the mask from frame 1. If the cell moves, shrinks, the mask should be made for each frame. 

  >> [avgfusioncellintensity_time] = CellIntensityInTime(fusion_movie, cellMask);
  >> [avgproteincellintensity_time] = CellIntensityInTime(protein_movie, cellMask);
5. CALCULATE THE AVERAGE (space and time) CELL AND BACKGROUND INTENSITY
  >> [averageCellIntensityWholeMovieFusion] = mean(avgfusioncellintensity_time);
  >> [averageCellIntensityWholeMovieProtein] = mean(avgproteincellintensity_time);
6.	FIND FUSION LOCATIONS: This will calculate a difference movie, then bandpass filter. Watch your movie first to have some idea of the number of events. Is it 1 or 100? This will help you set a threshold. You will need bpass, pkfind.
  >> thresholdFactor = 0.5; 
  %vary this threshold, I prefer to find non-fusion events instead of missing fusion events. In the code, this is multiplied by the max intensity, so it is a scaling factor and can accommodate different expression levels. 
  >> exposuretimems = 50; 
  %we stream at 50ms/frame and the fusion event finder code calculates a difference movie with a delta time of 0.5s. This will depend on how long events last and 0.5-1.25s works well for exosome secretion. The next command outputs XY position of potential fusion locations. 
  >> [fusion_locations] = FusionEventFinderV5(fusion_movie, cellMask, thresholdFactor, exposuretimems);
7.	OUTPUT MINI FUSION MOVIES and ALIGNED INTENSITY TRACES: Once you have locations, crop the raw movie file to make mini movie sequences (25x25 pixels) around each location. These will be aligned in time for frame 50 to be the fusion point. If fusion happens in a frame < 50, empty frames are added at the start. The cell mask is used to normalize the fusion channel: (circle – cell avg intensity) / (max – cell avg intensity). An annulus (local background) is used to normalize the protein channel, this is dF/dFmax:(circle – annulus) / (max – annulus).   %this calls MiniStackFusionDataProcessV9.m
This makes the following figures and outputs: 1) Show the onset of fusion point, the red dot is the zero point in time, 2) Aligned full traces for fusion, 3) aligned full traces for the protein channel. All aligned traces are written into excel (time, fusion, protein). 1 and 2 are shown below.
  >> [fusion_mini, protein_mini, data] = FusionMiniMaker8(fusion_movie, protein_movie, fusion_locations,'3242F', '3242P', averageCellIntensityWholeMovieFusion);

8.	SAVE WORKSPACE
9.	COMPILE ALL EXCEL INTO ONE EXCEL AND PLOT

 
