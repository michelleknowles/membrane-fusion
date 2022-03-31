function [cellMask, backgroundMask, averageCellIntensity] = CellMaskMaker3(fusionMovie, ThresholdFactor, filename)
%CALL: [cellMask, backgroundMask, avgcellintensity] = CellMaskMaker2(fusion_movie, 0.3, 'masks.tif');
%INPUT: fusion movie, called from a variable in your workspace
%       ThresholdFactor - this is a factor that dictates the intensity that
%       is used to draw the border.
%       filename: the name of the output file that is a movie 
%OUTPUT: a 2 x 256x256 array that are the cell mask and the background mask.
%Adapted from Anarkali's Fusion Process code. AW 11/8/21
%BB 3/1/22 to call a variable instead of reading in a movie.
%MKK 3/16/22 streamlining workflow and commenting. 

%look at the first image for making the mask
imageForCellMaskArray = fusionMovie(:,:,1);

%Obtain cell and background masks
[cellMask,backgroundMask] = CreateCellMask(imageForCellMaskArray, ThresholdFactor);

%APPLICATION OF CELL MASK
%Calculate relative expression level(RExLevel) = cell intensity in the first image minus background
imageForCellMaskArray = double(imageForCellMaskArray); 
%Multiply by cell mask so intensity value outside cell will become zero
appliedCellMaskSingleFrame=imageForCellMaskArray.*cellMask;

% Calculate Avg cell intensity in first image of the cell
averageCellMaskSingleFrame = mean(nonzeros(appliedCellMaskSingleFrame));
averageCellIntensity = averageCellMaskSingleFrame;

%%APPLICATION OF BACKGROUND MASK
% Calculate Avg background in first image of the cell using background mask
backgroundCellMaskSingleFrame = imageForCellMaskArray.*backgroundMask;
averageBackgroundSingleFrame = mean(nonzeros(backgroundCellMaskSingleFrame)); 
    
%DISPLAY MASK TO USER% 
contrastHigh1 = averageCellMaskSingleFrame * 2.5;
contrastLow = averageBackgroundSingleFrame * 0.5;
contrastHigh2 = averageCellMaskSingleFrame * 1.2;
    
figure;
imshow(appliedCellMaskSingleFrame, [contrastLow,contrastHigh1]);
  
figure;
imshow(backgroundCellMaskSingleFrame, [contrastLow,contrastHigh2]);
     
%%WRITE DATA/OUTPUTS
    path = pwd;
    outputPath = [path '\'];
            
%Display the mask 
appliedCellMaskSingleFrame = uint16(appliedCellMaskSingleFrame);
backgroundCellMaskSingleFrame = uint16(backgroundCellMaskSingleFrame);

%%Name the masks
fusionMovieMASKNameAppended = append(filename,'_CellMaskImage');
fusionMovieBACKGROUNDNameAppended = append(filename,'_BackgroundMaskImage');

% ARRAY DATA OUTPUT
% these are binary masks 256x256x1
cellMask = cellMask;
backgroundMask = double(backgroundMask);
backgroundMask = backgroundMask;
    
%Output array information about cell and background intensity.

end

function [cellmask,backgroundmask]=CreateCellMask(cellMaskImageArray, ThresholdFactor)
% 
% PURPOSE:  to create cell mask from a single image
% INPUT: an image, a threshold
% OUTPUT: Cellmask and backgroundmask. 
% AM built code with help from https://www.mathworks.com/help/
% MKK documented 3/16/22

    %This controls the intensity that gets cutout. The lower it is, the
    %more cell material you will pickup. If it is 1, you cut out almost
    %everything that isn't max intensity...
    intensityThresholdFactor = ThresholdFactor;
    
    %Taken from: https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
    %Find the edge, requires the Image Processing Toolbox
    [~, threshold] = edge(cellMaskImageArray, 'sobel');
    BWs = edge(cellMaskImageArray,'sobel',threshold * intensityThresholdFactor);
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    BWsdil = imdilate(BWs,[se90 se0]);
    BWdfill = imfill(BWsdil,8,'holes');
    BWnobord = imclearborder(BWdfill);
    bckcormsk=~BWnobord(:,:);
    cellmask = BWdfill-BWnobord;
    backgroundmask(:,:) = ~cellmask(:,:);
    seD = strel('diamond',1);
    BWfinal1 = imerode(cellmask,seD);
    cellperimeter = bwperim(BWfinal1);
    Segout = cellMaskImageArray;
    Segout(cellperimeter) = 255;
% 
end

    
    
    
    
    