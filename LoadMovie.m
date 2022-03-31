function movie_out = LoadMovie(movie)
%I want this script to look at the movie being loaded and determine if it
%is a .tif or a .stk and then load it appropriately. Convert to the
%appropriate format to work with.


    %possible extensions
    tif_extension = '.tif';
    stk_extension = '.stk';
    
    %getting data from the movie
    [filepath,name,ext] = fileparts(movie);
    movie_ext = ext;
    
    %loading the movie depending on extension
    if movie_ext == tif_extension
        LoadedMovie = TIFFStack(movie);
    elseif movie_ext == stk_extension
        LoadedMovie = stkread(movie, filepath);
    else
        disp('incompatible movie');
    end
    
    %convert movie to uint16
    LoadedMovie = LoadedMovie(:,:,:,1);
    
    movie_out = LoadedMovie;
end