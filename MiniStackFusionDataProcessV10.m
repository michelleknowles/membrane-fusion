function out = MiniStackFusionDataProcessV10(fusionMiniStack, proteinMiniStack, fusionMovieName, proteinMovieName, averageCellIntensityInTime, averageProteinIntensityInTime)
% PURPOSE: This takes a fusion and protein ministack movies paired together and does the following:
%  1) multiplies by a mask
%  2) measures the average intensity in a circle in the center. 
%  3) The intensity traces have the background subtracted: the average cell
% intensity in time is subtracted for fusion and the annulus is subtracted
% for the protein channel. These are different because the fusion event
% clearly secretes into the annulus area. 
%  4) The intensity traces are aligned in time so that the onset is 0s.
%  This is done by looking at the fusion event then aligning both the
%  fusion and protein. 
%  5) The onset point is plotted in one figure and the entire aligned trace
%  is plotted in another. 
% INPUT: 
% fusionMiniStack: a cropped 25x25 pixel movie of a potential fusion event. Note the 
% FusionMiniMaker code outputs a 4D variable, choose which spot you want to
% analyze BEFORE inputting (e.g. fusionMiniStack(:,:,:,x) where x is spot of
% interest)
% proteinMiniStackFull: same notes as fusionMiniStackFull
% fusionMovieName and proteinMovieName: choose output names for your files.
% averageCellIntensityInTime: This is the entire cell intensity measured
% after applying a mask.
% OUTPUT: an intensity file with the time (aligned to onset 0s), intensity of fusion, intensity of protein 
% NOTE: the normalization for fusion and protein are different so if you
% send in the fusion movie for both channels, you will get different but
% similar traces out. 
% CHANGES: Taken from V9 on 4/20/23 MKK Adapted to fix a concatentation
% error. Catch/try in FusionMiniMaker probably not needed now. 

%Create mask
[MaskC MaskA] = maskMaker3(7,11);
MaskC = double(MaskC);
MaskA = double(MaskA);
nFrames = size(fusionMiniStack,3);
MovieLength = nFrames;

%pad the movies with NaN frames and correspondingly shift the average
%cellintensityintime. 
preFusionFrames = 50;
NaNframes = NaN(25, 25, preFusionFrames);
NaNpaddedFusionmovie = cat(3, NaNframes, fusionMiniStack);
NaNpaddedProteinMovie = cat(3, NaNframes, proteinMiniStack);
fusionMiniStack = NaNpaddedFusionmovie;
proteinMiniStack = NaNpaddedProteinMovie;
nFrames = nFrames + preFusionFrames;
NaNpadIntensity = NaN(preFusionFrames,1);
averageCellIntensityInTime = cat(1,NaNpadIntensity,averageCellIntensityInTime);
averageProteinIntensityInTime = cat(1,NaNpadIntensity,averageProteinIntensityInTime); 

%Multiply both MiniStackMovies by the MiniStackMask
fusionMiniStackC = (double(fusionMiniStack).*MaskC);
proteinMiniStackC = (double(proteinMiniStack).*MaskC);
fusionMiniStackA = (double(fusionMiniStack).*MaskA);
proteinMiniStackA = (double(proteinMiniStack).*MaskA);

%Find the Average Intensity for each frame of the fusionMiniStack Movie
[fusionIntensityArrayFull] = AverageIntensityPerFrame(fusionMiniStackC, MaskC);
[fusionIntensityArrayFullBG] = AverageIntensityPerFrame(fusionMiniStackA, MaskA);
fusionIntensityArrayAvgBG = mean(fusionIntensityArrayFullBG,3,'omitnan');
fusionIntensityArrayAvgBG = mean(nonzeros(fusionIntensityArrayAvgBG),'omitnan');

%Find the frame with maximum Intensity and store it
[maxIntensityOfFusionMovieFrame] = FusionIntensityMaxDetector(fusionIntensityArrayFull);

%If the max intensity - 50 is less than 2 change it to 3
if (maxIntensityOfFusionMovieFrame - 50) < 3
    maxIntensityOfFusionMovieFrame = 53;
else
end
    

%normalize the fusionIntensityArrayFull
[normalizedFusionIntensityArrayFull] = FusionDataNormalize(fusionIntensityArrayFull, averageCellIntensityInTime,nFrames);

%Select the prefusion frames to do analysis on.
%MKK changed to start at preFusionFrames rather than 1 and it worked for some but not all traces
preFusionNormalized = normalizedFusionIntensityArrayFull(1:maxIntensityOfFusionMovieFrame);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Take the max intensity frame number and then select all frames leading up
%to that and separate. This is post normaization. Then run the ischanged
%function with the various thresholds and parameters on this prefusion
%leadup.

%This is where I will put the prefusionFrameDetector using difference in
%intensity to detect the spike.
[preFusionIndex] = PreFusionAlignment(preFusionNormalized, preFusionFrames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get total frames of the movie and set the post frames as maximum minus the
%fusion event frame.
miniStackMovieLengthFull = size(fusionMiniStack);
miniStackMovieLengthFull = miniStackMovieLengthFull(1,3);
%This is the whole movie length subtracted by the high intensity event
%frame.
postFusionFrames = miniStackMovieLengthFull - preFusionIndex;

%These are the selected frames from the full movies.
[fusionMiniStackMovieSelect, NewPreFrames] = FusionMiniMovieSelect(fusionMiniStackC, preFusionIndex, preFusionFrames, postFusionFrames);
proteinMiniStackMovieSelect = FusionMiniMovieSelect(proteinMiniStackC, preFusionIndex, preFusionFrames, postFusionFrames);

fusionMiniStackMovieSelectBG = FusionMiniMovieSelect(fusionMiniStackA, preFusionIndex, preFusionFrames, postFusionFrames);
proteinMiniStackMovieSelectBG = FusionMiniMovieSelect(proteinMiniStackA, preFusionIndex, preFusionFrames, postFusionFrames);

%Get intensity of these frames %This is the average intensity of the
%frame
[fusionIntensityArraySelect] = AverageIntensityPerFrame(fusionMiniStackMovieSelect, MaskC);
[proteinIntensityArraySelect] = AverageIntensityPerFrame(proteinMiniStackMovieSelect, MaskC);

[fusionIntensityArraySelectBG] = AverageIntensityPerFrame(fusionMiniStackMovieSelectBG, MaskA);
[proteinIntensityArraySelectBG] = AverageIntensityPerFrame(proteinMiniStackMovieSelectBG, MaskA);


%THIS IS AN ATTEMPT TO CORRECT FOR THE FAST PHOTOBLEACHING OF RED. May not
%be needed.
%[normalizedDataArrayFusion] = fusionIntensityArraySelect-fusionIntensityArraySelectBG;
%[normalizedDataArrayProtein] = proteinIntensityArraySelect-proteinIntensityArraySelectBG;

%Data normalization, be careful with NaNs here
%[normalizedDataArrayFusion] = FusionDataNormalize(fusionIntensityArraySelect, averageCellIntensityInTime,nFrames);
%[normalizedDataArrayProtein] = FusionDataNormalize(proteinIntensityArraySelect, averageProteinIntensityInTime,nFrames);

%THIS CORRECTS FOR THE FACT THAT the normalized intensity prior to fusion is too high. We need
%to subtract a local background (A Annulus)rather than the entire cell, but the annulus increases post
%fusion due to release. Therefore this uses the annulus intensity prefusion only. 

%MaskA first frame for both protein and fusion.
%maskAfusionFirstFrame = fusionIntensityArraySelectBG(1,1);
%maskAproteinFirstFrame = proteinIntensityArraySelectBG(1,1);
%calc the average of the annulus prefusion (just 5 frames) This is at time
%-50 to -46 frames before max. Matches what Brody does. 
maskAfusionFirstFrame = mean(fusionIntensityArraySelectBG(1:5,1))
maskAproteinFirstFrame = mean(proteinIntensityArraySelectBG(1:5,1))
%%%%%%%%%
%make array that is size of full movie typically 550 frames but make it
%versatile and fill with those numbers then feed into the normalization

sz = [nFrames 1];

emptyArrayN = ones(sz);
%this makes arrays filled with one number - the average annulus pre-fusion
fusionAarray = emptyArrayN * maskAfusionFirstFrame;
proteinAarray = emptyArrayN * maskAproteinFirstFrame;

[normalizedDataArrayFusion] = FusionDataNormalize(fusionIntensityArraySelect, fusionAarray,nFrames);
[normalizedDataArrayProtein] = FusionDataNormalize(proteinIntensityArraySelect, proteinAarray,nFrames);

%[normalizedDataArrayProtein] = ProteinDataNormalize(proteinIntensityArraySelect, proteinIntensityArraySelectBG);

%%% WRITE FILES %%%
path = pwd;
outputPath = [path '\'];
            
fusionStkFileName = [fusionMovieName,' Mask','.stk'];
proteinStkFileName = [proteinMovieName,' Mask','.stk'];

fusionMiniStackMovieSelect = uint16(fusionMiniStackMovieSelect);
nFramesMask = size(fusionMiniStackMovieSelect,3)
proteinMiniStackMovieSelect = uint16(proteinMiniStackMovieSelect);

%stkwrite(DATA,NAME,outputPath);
if  nFramesMask < MovieLength
    stkwrite(fusionMiniStackMovieSelect,fusionStkFileName,outputPath);
    stkwrite(proteinMiniStackMovieSelect,proteinStkFileName,outputPath);
else
    stkwrite(fusionMiniStackMovieSelect(:,:,1:MovieLength),fusionStkFileName,outputPath);
    stkwrite(proteinMiniStackMovieSelect(:,:,1:MovieLength),proteinStkFileName,outputPath);
end

%Output the miniStack.stk files for each.

%Output .xls for each I could put them in the same file, but I will not for
%now.

%I could append some other set of numbers that works with the pre post
%fusion events so that everything is already labeled based on frames. Or
%even ms time Append so that both columns are in the same sheet.

%Creates a column of numbers with the peak intensity of the fusion event as
%0. 
[fusionArrayFrames] = FrameArrayCreator(NewPreFrames, postFusionFrames);
if length(fusionArrayFrames) > 500 
    fusionArrayFrames = fusionArrayFrames(1:500);
end
%xlsWrite with 3 columns I can use Array Size to work on this too...
%Combine data into a single array and write as a range?

%[~,name]=fileparts(fusionMovieName);
%xlsWriteNameFusionRemoveExtension = name;
%clear name;

%xlsWriteNameFusion = [xlsWriteNameFusionRemoveExtension,' IntensityTrace'];


%First Column is the timeArray, Second is the Fusion Data, Third is
%the Protein Channel.
%CHANGED IN V10: This needed a check for movies that last 500 frames, the
%fusionArrayFrames is 1 too long to concatenate. 
out = [fusionArrayFrames normalizedDataArrayFusion normalizedDataArrayProtein];

%xlswrite(xlsWriteNameFusion, fusionArrayFrames,'Sheet1','A');
%xlswrite(xlsWriteNameFusion, normalizedDataArrayFusion,'Sheet1','B');
%xlswrite(xlsWriteNameFusion, normalizedDataArrayProtein,'Sheet1','C');

%I am going to do some checks here to make sure things are working
%correctly.
%What is the max intensity frameIndex?
maxIntensityOfFusionMovieFrame
%What is the preFusionIndex?
preFusionIndex
%What is the postFusionSize? This includes the fusion event
postFusionFrames
%Do things add up correctly and is it the same as our output?
preFusionFrames + postFusionFrames %This is likely missing a single frame, maxIntensityOfFusionMovieFrame
%What is the size of the output? 
length(out)

preFusionIndex + postFusionFrames

%What is the new maxIntensity after the prefusion cropping has happened?
preFusionFrames + (maxIntensityOfFusionMovieFrame - preFusionIndex)
%What is the prefusionIndex after the cropping?
    %This should be 50

%%PLOTTING%%
%PLOT 1: to show onset of fusion as a red circle
figure(1)
plot(preFusionNormalized);
xlabel('time(frames)');
ylabel('intensity');
hold on
%identify the intensity spike onset and plot it
newZero = preFusionNormalized(preFusionIndex,1);
newZero = [preFusionIndex, newZero];
plot(newZero(:,1), newZero(:,2), 'ro');
%PLOT 2: show entire aligned fusion trace 
figure(2)
xlabel('time(frames)');
ylabel('intensity');
hold on
plot(fusionArrayFrames, normalizedDataArrayFusion);
%PLOT 3: show entire aligned protein trace
figure(3)
xlabel('time(frames)');
ylabel('intensity');
hold on
plot(fusionArrayFrames, normalizedDataArrayProtein);

end

%%%%% FUNCTIONS USED %%%%%
function [averageIntensityPerFrameTotal] = AverageIntensityPerFrame(miniStackMovieMaskMultiply, maskImage)

%PURPOSE: This converts the mini movie into an intensity in time for a
%circular region in the center of the cropped mini movie. 
%INPUT: Mini movie (25x25 pixels, nframes) and a mask of the desired size
%OUTPUT: Intensity of the masked region in time. 
%
%AW 11/4/21

averageIntensityPerFrameTotal = [];
maskSize = size(miniStackMovieMaskMultiply);
frames = maskSize(3);
maskSize = sum(maskImage(:));

for i = 1:frames
    
    singleFrame = miniStackMovieMaskMultiply(:,:,i);
    frameIntensity = sum(singleFrame(:));
    averageIntensitySingleFrame = (frameIntensity/maskSize);
    
    averageIntensityPerFrameTotal = [averageIntensityPerFrameTotal;averageIntensitySingleFrame];
end

averageIntensityPerFrameTotal = averageIntensityPerFrameTotal; %output
end

function [maxFusionFrameNumber] = FusionIntensityMaxDetector(intensityArray)

%PURPOSE: To find the maxiumum intensity in the intensity array.
%HOW: Takes an array of the average intensity of a ministack movie, from the
%function AverageIntensityPerFrame and finds the maximum value and reports
%back the frame# or the array coordinates that the max intensity occurs at.
%INPUT: an intensity trace in time.
%OUTPUT: a location in the array where the maximum occurs.

[~, ind] = max(intensityArray);
frameNumber = ind;
maxFusionFrameNumber = frameNumber;
end

function [selectFusionMiniStack, NewPreFrames] = FusionMiniMovieSelect(miniStackMovie, maxFusionFrameIntensity, preFusionFrames, postFusionFrames)
%PURPOSE: To make a minimovie that is a subset of frames from a larger
%movie where the fusion event position is identified and will be assigned
%to a specific frame number (frame# = prefusionFrames)
%INPUT: a single mini movie of a fusion event
%preFusionFrames are the number of frames leading up to fusion,
%postFusionFrames are the frames that are post fusion.
%OUTPUT: Selects this section of frames and makes a new miniStackMovie out of them.
%Writes mini movie as a stack. 
%NOTE: If an event shows up on frame 10 and the number of prefusion frames
%is larger, the event is still kept and the beginning frames are NaN
%filled. 
%AW 11/4/21

%makes a series of frames filled with NaNs; preFusionFrames (50) is a constant for all movies 
%NaNframes = NaN(25, 25, preFusionFrames);
%concatenate NaN frames to ensure there are at least 50 frames prior to
%fusion. (1050 frames now with 50 NaN frames at the start of the movie.
%This means the maxfusionframe is 50 higher.
%NaNpaddedmovie = cat(3, NaNframes, miniStackMovie);

NewPreFrames = preFusionFrames;
post = (maxFusionFrameIntensity + postFusionFrames);
pre = (maxFusionFrameIntensity - preFusionFrames);
selectFusionMiniStack = miniStackMovie(:,:,pre:post);

end

function [fusionArrayFrames] = FrameArrayCreator(preFusionFrames, postFusionFrames)

%This is going to take in some fusion event frame numbers and return  an
%array that can be appended to the xlswrite file so that there is an x axis
%for everything.
%AW 2021

negativeFrames = preFusionFrames * -1;

positiveFrames = postFusionFrames;

totalFrames = [negativeFrames:postFusionFrames];
totalFramesTransposed = transpose(totalFrames);

fusionArrayFrames = totalFramesTransposed;

end

function [normalizedDataArray] = FusionDataNormalize(miniMovieIntensityArraySelect, averageCellIntensityWholeMovie,nFrames)
%This takes in the array of data aligned with time and a float variable
%that is the average whole cell intensity and then normalizes with the
%below equation
% (circleArray - <cellIntensity from mask>) / (max(circle) - <cellIntensity from
% mask>)


%get max(circle)
maxIntensity = max(miniMovieIntensityArraySelect);
%maxIntensity = repmat(maxIntensity,1000,1);
size(maxIntensity);

nFrameSelect = nFrames+1-size(miniMovieIntensityArraySelect(:),1);
%miniMovieIntensityArraySelect=miniMovieIntensityArraySelect(nFrameSelect:nFrames);
averageCellIntensityWholeMovie = averageCellIntensityWholeMovie(nFrameSelect:nFrames);

size(averageCellIntensityWholeMovie);

size(miniMovieIntensityArraySelect);

%watch NaNs here!
diff1 = (miniMovieIntensityArraySelect - averageCellIntensityWholeMovie);
diff2 = (maxIntensity - averageCellIntensityWholeMovie);
normalizedDataArray = diff1./diff2;
normalizedDataArray = nonzeros(normalizedDataArray);

end

function [normalizedDataArray] = ProteinDataNormalize(miniMovieIntensityArraySelect, miniMovieIntensityBG)
%This takes in the array of data aligned with time and a float variable
%that is the average whole cell intensity and then normalizes with the
%below equation
% (circleArray - <cellIntensity from mask>) / (max(circle) - <cellIntensity from
% mask>)


%get max(circle)
maxIntensity = max(miniMovieIntensityArraySelect);

miniMovieIntensityBGmean = mean(miniMovieIntensityBG,3);
miniMovieIntensityBGmean = mean(nonzeros(miniMovieIntensityBGmean));



%
normalizedDataArray = (miniMovieIntensityArraySelect - miniMovieIntensityBG)./(maxIntensity - miniMovieIntensityBG);



end

function [preFusionFrameIndex] = PreFusionAlignment(intensityArray, preFusionFrames)
%figure out how to align to the base value, when something starts to change
%drastically, mark that as the 0 point... Should be pretty simple.

%This will take in the raw data and then give back a frame number? This is
%going to replace the smaxFusionFrameNumber that is fed into the
%selectFusionMiniStack



miniStackMovieLength = length(intensityArray);
%Why do we do this? The added NaNs will break the mean and STD calculations
%because they depend on indexing.
intensityArray = intensityArray(preFusionFrames:miniStackMovieLength);%preFusionFrames is not an index value 
miniStackMovieLength = length(intensityArray);


%%% This part is okay and I trust
quarterLengthIndex = round(miniStackMovieLength/4);
halfLengthIndex = round(miniStackMovieLength/2);
threeQuarterLengthIndex = round((miniStackMovieLength/4)*3);

miniStackQuarterMean = mean(intensityArray(1:quarterLengthIndex,1),'omitnan');
quarterSTD = std(intensityArray(1:quarterLengthIndex,1),'omitnan');
miniStackHalfMean = mean(intensityArray(1:halfLengthIndex,1),'omitnan');
halfSTD = std(intensityArray(1:halfLengthIndex,1),'omitnan');

miniThreeQuarterLengthMean = mean(intensityArray(halfLengthIndex:threeQuarterLengthIndex,1),'omitnan');
miniThreeQuarterLengthSTD = std(intensityArray(halfLengthIndex:threeQuarterLengthIndex,1),'omitnan');
%Get length of the movie.
%Get the initial mean.
%see how linear the initial half of the data points are.
%walk backwards through the data and test if the data point falls within
%the mean. 
%Once it does, mark that index.
%Acquire the mean and the STD of the data.
%%%


%Walk backwards through the data and mark when the data point falls within
%the expected mean/STD Only look at the last quarter of the data?
for i = miniStackMovieLength:-1:threeQuarterLengthIndex
    currentIntensity = intensityArray(i);
    %gather data of active index 
    
    rangeLower = (miniThreeQuarterLengthMean - miniThreeQuarterLengthSTD*3);
    rangeUpper = (miniThreeQuarterLengthMean + miniThreeQuarterLengthSTD*1);
    
    if currentIntensity > rangeLower && currentIntensity < rangeUpper
    %compare datapoint with the mean stuff.....
    %mean+STD
    
        preFusionFrameIndex = i + preFusionFrames - 1 %why do I have to subtract one?? Because you removed the NaN values not based on an Index value but on the number of NaNs we wanted to add.
    return;
        %set the index out variable and break out of the loop
    else
        %keep going through the loop
    end

end
%Possibly walk backwards from this and check that the point is within the
%mean of the first few points, some degree...

%This threshold changes where the 0 peak is found. At 1 it is as if it is
%not there. 
%detectChange = ischange(intensityArray, 'linear', 'Threshold', 4);

%possibly a better way to do this...
%firstChange = find(detectChange,1,'last');
end


