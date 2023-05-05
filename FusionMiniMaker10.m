function [fusionOut, proteinOut, dataOut] = FusionMiniMaker10(fusionChannelMovie, proteinChannelMovie, coords, fusionMovieName, proteinMovieName,averageCellIntensityWholeMovieFusion,averageCellIntensityWholeMovieProtein)

% PURPOSE: To crop out 25x25 pixel mini movies from the fusion sites, measure the intensity within a circle in the center, then find the onset of fusion.
% CALL: [fusion_mini, protein_mini, data] = FusionMiniMaker3(r,g,pk_event,'1356F','1356P',avgcellintensity_time);
% INPUT: fusionChannelMovie: fusion movie already in workspace
%   proteinChannelMovie: protein movie already in workspace. Use the fusion movie again if you do not have a second channel. 
%   coord: peaks found using FusionEventFinderV5
%   fusionMovieName: name to begin fusion movie stacks ('_x' will be appended, where x is the spot #, so specify fusion here (e.g. I use '####F', so an outputted filename would be '0000F_1' for location # 1 in cell ####)
%   proteinMovieName: same as fusion, but name for your protein channel
%   (e.g. '####P')
% OUTPUT: fusionOut outputs a 4D matrix in your variables of your fusion
% channel. The format of this is (25,25, nframes, nfusionevents)
% proteinOut likewise outputs a 4D matrix of your protein channel
% dataout = output of the aligned intensties in time for each location of a
% potential fusion event. Note these are also written to an excel file. 
% UPDATES:
% Adapted from FusionMiniMaker4 (AW) by BB on 3/1/22 to call variables from
% your workspace and output variables alongside written ministacks
% MKK 3/17/22 updated code for better commenting and removed unused lines.

%Make 4D array of ministacks, note the first movie is blank
fusionMiniMovie4D = ministk_movie(fusionChannelMovie,coords,25);
proteinMiniMovie4D = ministk_movie(proteinChannelMovie,coords,25);

%Get size of 4D array.
Size4D = size(fusionMiniMovie4D);
nframes = Size4D(1,3);
Size4D = Size4D(1,4);
%Might need to do some debugging here based on how large coord is.

fusionstk = zeros(25,25,nframes,Size4D);
proteinstk = zeros(25,25,nframes,Size4D);

for i = 2:Size4D %why does this start at 2?
    
    activeFusionMiniStack = fusionMiniMovie4D(:,:,:,i);
    activeProteinMiniStack = proteinMiniMovie4D(:,:,:,i);

    fusionstk(:,:,:,i) = activeFusionMiniStack;
    proteinstk(:,:,:,i) = activeProteinMiniStack;
    
    currentLoopAsString = num2str(i - 1);
    SpotNumber = currentLoopAsString
     
    fusionMovieNameAppended = append(fusionMovieName,'_',currentLoopAsString);
    proteinMovieNameAppended = append(proteinMovieName,'_',currentLoopAsString);
    
    %%%% WRITE FILES %%%%
    path = pwd;
    outputPath = [path '\'];
            
    fusionStkFileName = [fusionMovieNameAppended,'.stk'];
    proteinStkFileName = [proteinMovieNameAppended,'.stk'];
    %CALC INTENSITIES!!
    %This is where all the business happens: the mini movie is sent in to be multiplied by a mask (a central circle) and then outputs intensties for both channels. 
    %Each channel is normalized differently: fusion is normalized by the average cell intenstity in time and the protein is normalized by a local background annulus to give dF/S. 
    try
       datatemp = MiniStackFusionDataProcessV10(activeFusionMiniStack,activeProteinMiniStack,fusionStkFileName,proteinStkFileName,averageCellIntensityWholeMovieFusion,averageCellIntensityWholeMovieProtein);
    catch
        datatemp = [0 0 0];
        warning('Problem using function.  Assigning a value of 0.');
    end
    
    sizei = size(datatemp,1);
    sizeNaN = sizei+1;
    datacompile(1:sizei,:,i-1) = datatemp;
    datacompile(sizeNaN:nframes,:,i-1) = NaN;
    
    %stkwrite(DATA,NAME,outputPath);
    stkwrite(activeFusionMiniStack,fusionStkFileName,outputPath);
    stkwrite(activeProteinMiniStack,proteinStkFileName,outputPath);

    clear fusionMovieNameAppended;
    clear proteinMovieNameAppended;
    clear datatemp;
    
end

fusionOut = fusionstk(:,:,:,2:Size4D);
proteinOut = proteinstk(:,:,:,2:Size4D);
dataOut = datacompile(1:nframes,:,:);
writematrix(dataOut, 'alldecays.xls','Sheet','decays');
end

%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
function out=ministk_movie(im,rgn,sz)
% out=ministk(im,rgn,sz)
% 
% PURPOSE:  to cut out small regions of an image based on spots found in
% the image or a corresponding image of a different color. This is used to
% measure colocalization based on the work of Knowles and Barg in two PNAS
% 2011 papers. Spots are found using the work of tracking routines
% available on Eric Weeks' website (Emory University) and made into Matlab
% by Eric Dufrense.
% INPUT:
% im: movie stack to cut from
% rgn: spots (x,y) about which regions should be cut
% sz: size of cut out (a square of sz by sz pixels)
% sepdist: is the minimum separation distance between two spots. If two
% spots are within this distance of one another, neither are counted.
%
% OUTPUT:  a sz x sz x N_frames array for each spot
%           
% CREATED: Michelle Knowles May 2012
% EDITED: to crop movies rather than single images by Michelle Knowles and
% Aubrie Blevins December 2017, note that X and Y from the peak finder seem
% to be switched when cropping in this way.

%if sz/2 == floor(sz/2)
%warning('sz must be even so that the spots can be centered on a pixel: 1 pixel added');
%sz = sz+1
%end

%scott's code for a mask
pix=(sz+1)/2;
dimy = length(im(1,:,1));
dimx = length(im(:,1,1));
nrgn = length(rgn(:,1));
nframes = length(im(1,1,:));

%create a blank image array that you can fill
msk=zeros([sz,sz,nframes]);
%loop through all regions that locate spots in an image
for i =1:nrgn;
    x = rgn(i,1);
    y = rgn(i,2);
    %don't include regions within pix distance from the edge of the image.
    if ((x>pix) & ((x+pix)<dimy) & (y>pix) & ((y+pix)<dimx))
          cutout=im(y-pix+1:y+pix-1, x-pix+1:x+pix-1,:);
          msk = cat(4, msk, cutout);
    end
end
%out=msk(:,:,:,2:nrgn);
out=msk;
end


