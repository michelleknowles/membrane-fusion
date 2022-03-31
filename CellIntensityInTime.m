function [avgcellintensity_time] = CellIntensityInTime(cell_movie, cellMask);
%CALL: [avgcellintensity_time] = CellIntensityInTime(fusion_movie, cellMask);
%INPUT: a movie of a cell, typically of type uint16, and the cell mask that
%was made using CellMaskMaker3. 
%OUTPUT: the average cell intensity in time. Only the part of the image
%within the cell mask is included in the average. 

b=double(cell_movie);
maskedcell=b.*cellMask;
nframes = length(b(1,1,:));
avgcellintensity_time=zeros(nframes,1);

for i = 1:nframes
    avgcellintensity_time(i) = mean(nonzeros(maskedcell(:,:,i)));
end