function [fusionCoordinatesOut] = FusionEventFinderV5(fusionMovie,cellMask, thresholdFactor, timePerFrameMS)

% CALL:[fusionCoordinatesOut] = FusionEventFinderV4(fusion_movie,cellMask,
% 20, 50);
% INPUT: fusionMovie: fusion channel, called into the workspace
%     cellMask: cell mask found by CellMaskMaker3
%     thresholdFactor: threshold factor to find peaks. It is multiplied by the maximum to get a threshold in intensity units. Increase to make less sensitive
%     and find brighter events only.
%     timePerFrameMS is the exposure time used when taking data. This is
%     used for calculating a difference movie. 
% OUTPUT: XY coordinates of possible fusion events
% USES: bpass and pkfnd from https://site.physics.georgetown.edu/matlab/
% UPDATES: 
% BB  3/1/22 Adapted from FusionEventFinder (AW) by to use a variable
% instead of reading in a movie
% AW 3/1/22 Changed Bpass filter order with calcDifferenceMovie.
% instead of calling a movie
% AW 3/8/22 dt based on timePerFrame 
% MKK 3/17/22 commenting 

dtTotalTime = 500; %standard ms time for the difference movie calculation
%dt depends on the image exposure time. This is rounded.
dt = ceil(dtTotalTime/timePerFrameMS); %For 50ms this should be dt=10 .

% Generate difference movie and then the maximum projection image of the
% difference movie
[~,maxProjectionDifferenceArray] = CalcDifferenceMovie(fusionMovie, dt);

% Bandpass filter the difference movie
filter = 7;
maxProjectionDifferenceArrayBPASS = bpass(maxProjectionDifferenceArray,1,filter);

% Multiply Max projection and cell mask to eliminate background
maxProjectionDifferenceArrayMasked = maxProjectionDifferenceArrayBPASS.*cellMask;

%Amplify signal
%Intensity values within cell perimeter are >0
cell_Intensity_values = (maxProjectionDifferenceArrayMasked)>0;

%Increase brightness of the bright spots for enhanced contrast and
%easier peak selection
multipliedMaxProjection = maxProjectionDifferenceArrayMasked.*cell_Intensity_values;

% enanced intensity values located within cell perimeter
enhancedMaxProjection = (nonzeros(multipliedMaxProjection));

contrastHigh = max(max(maxProjectionDifferenceArrayMasked));
contrastLow = 0;
    
%Find threshold to detect fusion event locations using pkfnd
%Set threshold using average brigntness of the cell
pkfndThresh = max(max(max(maxProjectionDifferenceArrayMasked)));
pkfndThresh = pkfndThresh * thresholdFactor;

%find XY location of bright spots
pkfndSize = 7; 
fusionCoordinates = pkfnd(multipliedMaxProjection, pkfndThresh, pkfndSize);

%plot what you found
figure;
imshow(maxProjectionDifferenceArrayBPASS, [contrastLow,contrastHigh]);
hold on
plot(fusionCoordinates(:,1), fusionCoordinates(:,2),'ro');

fusionCoordinatesOut = fusionCoordinates;
end
 

function [mov_out,maxproj] = CalcDifferenceMovie( mov, dt )
%PURPOSE: to calculate the difference in intensity from one frame to the
%next, output is a movie of the difference images.
%   CALL: [diffMovie, maxProj] = CalcDifferenceMovie(a, 20);
%   INPUT: a is a movie.  The format is usually (256,256, nframes)
%   dt is the time lag in frames over which the difference movie is calculated 
%   OUTPUT: this calculates a difference movie and then a maximum
%   projection of that movie
%   Steps - calculate difference images (frame t+1) - (frame t).
%  MKK Dec 2 2019
%
nframes = max(size(mov(1,1,:)));
x = max(size(mov(:,1,1)));
y = max(size(mov(1,:,1)));

mov_out = zeros(x, y, nframes-dt);
for i= 1:nframes-dt
    im = mov(:,:,i);
    im2 = mov(:,:,i+dt);

    % take the difference image from frame t+1 and frame t to determine new
    % pixels with values. Bright spots are where signal moved TO. 
    delta = im2 +1000 - im;
    %store the image in a movie file
    mov_out(:,:,i) = delta;

end

maxproj = max(mov_out, [], 3);
end
