%% This matlab code does all the plotting needed after icecap model runs (and more)
%% you just need to change study_directory accordingly to the right model run number folder
%% you might also need to change the DEM source files if running on different DEM
%% make sure that pathways to shapefiles are correct and havent't been modified 
%% for volume and area statistics, make sure the cell area is correct according to x and y cell sizes for a given DEM. This is calculated by projecting the DEM to the right UTM zone and reading x and y cell sizes in ArcGIS
%% for velocity maps and statistics: make sure the threshold value (used to get rid of high value anomalies) isn't too low: otherwise you will loose some output data.
%% for .tif rasters creation codes: make sure the corner coordinates are right. DEM should be geographic projection (UTM WGS84) and cropped to square (length = width) polygon for this to work well. 
%% you can check corner coodinates using this command and the source DEM .tif raster file:          info = geotiffinfo('namefile.tif');
%%                                                                                                  [x,y]=pixcenters(info); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


workenvironment = 'I:\tanc\tanc_field\Model_runs\Cordillera Blanca\Cojup_valley\Holocene_simulations\Model run 37\'; %% USER INPUT FOR WORKING DIRECTORY !! this needs to be changed for each model run 

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ice_elevation plots for data visualisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

study_directory = workenvironment;
files = dir(strcat(study_directory, 'thickness_*y.txt'));
domain = 'I:\tanc\tanc_field\domains\ALOS world 3d DEM\Cordillera_blanca\Cojup_glacier\WGS84\30m\no_modern_glacier\Lake_bump_trimmed\lake_removed_manually\larger_dem_for_holocene_simulations\LIA_moraine_removed\domain.dat';

for k=1 : length(files)

  %read the thickness file
  h0 = dlmread(strcat(study_directory, files(k).name));                 % read thickness files
  topo = dlmread(domain);                                               % read dem
  icesurface = topo + h0;                                               % add dem and ice thickness
  filename = strrep(files(k).name, 'thickness_', 'ice_elevation_');     % Replace in the filename 'thickness_' by 'ice_elevation_' 
  filename = strsplit(filename,'.txt');                                 % remove .txt extension from filename
  filenamepng = strcat(filename{1}, '.png');
  contourf(icesurface,60);                                              % make plot
  colorbar                                                              % add colorbar
  saveas(gcf,filenamepng)                                               % save plot as png
  clf                                                                   % clear figure

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ice_thickness plots for data visualisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


study_directory = workenvironment;
files = dir(strcat(study_directory, 'thickness_*y.txt'));

for k=1 : length(files)

  %read the thickness file
  h0 = dlmread(strcat(study_directory, files(k).name));                 % read thickness files
  filename = strrep(files(k).name, 'thickness_', 'ice_thickness_');     % Replace in the filename 'thickness_' by 'ice_elevation_' 
  filename = strsplit(filename,'.txt');                                 % remove .txt extension from filename
  filenamepng = strcat(filename{1}, '.png');
  contourf(h0,10);                                                      % make plot
  colorbar                                                              % add colorbar
  saveas(gcf,filenamepng)                                               % save plot as png
  clf                                                                   % clear figure

end


%%%%%%%%%%%%%%%%%%%%%%%%%% tif file ice elevation for opening in ArcGIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 h220 = dlmread('thickness_380y.txt');                                 % alter the Yr numbers in textfile title if your simulation reached steady-state earlier and was stopped
 topo = dlmread('domain.dat');
 icesurface = topo + h220;
 contourf(icesurface,60);
 colorbar

% Get geo referenced 
R = georasterref('RasterSize',size(topo),'LatitudeLimits',[-9.424437231036740,-9.358275393256740],'LongitudeLimits',[-77.403295559681750,-77.337133721901750]);

% write to tiff file 
tiffile = 'ice_elevation_360yr.tif' ;                                 % do not alter this title if you did alter the textfile input title. This is default used in the AAR python code... 
geotiffwrite(tiffile,icesurface,R)  

% you can now load the .tif file into ArcGIS if you want 
% reproject to UTM if you want to make measurements (otherwise in degrees and not metric)

%%%%%%%%%%%%%%%%%%%%%%%%%%% tif file ice thickness for opening in ArcGIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read thickness output files in folder 

h220 = dlmread('thickness_380y.txt');                               % alter the Yr numbers in textfile title if your simulation reached steady-state earlier and was stopped

% Get geo referenced 
R = georasterref('RasterSize',size(topo),'LatitudeLimits',[-9.424437231036740,-9.358275393256740],'LongitudeLimits',[-77.403295559681750,-77.337133721901750]);

% write to tiff file 
tiffile = 'glacier_thickness_360yr.tif' ;                           % do not alter this title if you did alter the textfile input title. This is default used in the AAR python code... 
geotiffwrite(tiffile,h220,R)  

% you can now load the .tif file into ArcGIS if you want 
% reproject to UTM if you want to make measurements (otherwise in degrees and not metric) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cajup glacier only thickness plots + volume/area/thickness statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %This code extracts the ice thickness output data to include only the data
%within a polygon imported as a shapefile, which was previously drawn in
%ArcGIS to highlight a specific area: currently the El Loro glacier catchment
%All that needs to be set is the location and name of the Domain (tiffile),
%the location and name of the polygon (shapefile), and the location and
%name of the thickness files (files) and its correct folder (study
%directory). The code then computes the volume for each thickness file
%present in the folder and also produces a PNG image for each of them.

shapefile = 'I:\tanc\tanc_field\Model_runs\Cordillera Blanca\Cojup_valley\Holocene_simulations\Shapefiles\glacier_polygon\Cojup_glacier_polygon.shp';
tiffile = 'I:\tanc\tanc_field\domains\ALOS world 3d DEM\Cordillera_blanca\Cojup_glacier\WGS84\30m\no_modern_glacier\Lake_bump_trimmed\lake_removed_manually\larger_dem_for_holocene_simulations\LIA_moraine_removed\Moraine_removed_DEM.tif';

polygon = shaperead(shapefile);
% When read, polygon.X, polygon.Y have trailing NaN, we need to get rid of
% it
rx = polygon.X(1: end-1); % take all array elements except the last one (NaN)
% same for Y:
ry = polygon.Y(1: end-1);

% Read tif file
geo = geotiffinfo(tiffile);
% Get the pixel coordinates (convert pixels to X,Y)
[x,y] = pixcenters(geo);                                                  % x and y are vectors
% Convert vectors to matrices
[X, Y] = meshgrid(x, y);
% Flip upside down the Y matrix because of the tif file...
Y = flipud(Y);

% Generate the mask: confront the polygon with the tiffile to get what is
% inside the polygon and what is outside
mask = inpolygon(X,Y, rx, ry);                                            % mask is a logical (1,0) matrix
mask = double(mask);                                                      % convert logical to double so as to multiply it with the elevation

% Read all thickness files
study_directory = workenvironment;
files = dir(strcat(study_directory, 'thickness_*y.txt'));
volume = [];

for k=1:length(files)
    % Read the thickness file
    h0 = dlmread(strcat(study_directory, files(k).name));

    %masked h0:
    res = times(h0, mask);                                                    %multiply each cell by its equivalent cell in other matrix
    contourf(res, 15);
    colorbar;
    volume_cojup(k) = 683.172669009681000* sum(sum(res));                     % 683...... value will change depending on DEM spatial resolution and thus cell area: creates an array indicating ice volume for each time step saved
    icearea_cojup(k) = sum(sum(res>0))*683.172669009681000;                   % will create another array indicating ice areal cover for each time step saved 
    max_thickness_per_file(k) = max(res(:));                                  % will create another array indicating maximum ice thickness over cropped glacier for each time step saved
    mean_thickness_per_file(k) = mean(res, 'all','omitnan');                  %calculates mean thickness across the glacier, thus for the entire polygon-cut matrix ('all'), and by omitting the NaN values ('omitnan').
    filename = strrep(files(k).name, 'thickness_', 'Cojup_Glacier_only_');    % Replace in the filename 'thickness_' by 'Cojup_glacier_only_' (nothing)
    filename = strsplit(filename,'.txt');                                     % remove .txt extension from filename
    title(filename, 'Interpreter', 'none');                                   % Interpreter: none means that the _ will not be interpreted as a subscript
    filenamepng = strcat(filename{1}, '.png');
    saveas(gcf, filenamepng)                                                  % export current figure (gca) into a .png file rather than a .fig file
    clf % clear figure ?
end

max_volume = (max(volume_cojup(:)))*10^-9;                                  %output glacier volume in km3
max_area = (max(icearea_cojup(:)))*10^-6;                                   %output glacier area in km2
max_thickness = max(max_thickness_per_file(:));                             %this calculates maximum ice_thickness for steady state glacier
mean_thickness = max(mean_thickness_per_file(:));                           %this calculates mean glacier ice_thickness for steady state glacier.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cojup glacier surface gradient code: to generate glacier elevation profile graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Open shapefile:

shapefile = 'I:\tanc\tanc_field\Model_runs\Cordillera Blanca\Cojup_valley\Holocene_simulations\Shapefiles\centreline\Cojup_centreline.shp';
tiffile = 'I:\tanc\tanc_field\domains\ALOS world 3d DEM\Cordillera_blanca\Cojup_glacier\WGS84\30m\no_modern_glacier\Lake_bump_trimmed\lake_removed_manually\larger_dem_for_holocene_simulations\LIA_moraine_removed\Moraine_removed_DEM.tif';

flowline = shaperead(shapefile);

%get rid of NaN

rx = flowline.X(1: end-1);
ry = flowline.Y(1: end-1);

% Read tif file
geo = geotiffinfo(tiffile);
% Get the pixel coordinates (convert pixels to X,Y)
[x,y] = pixcenters(geo);                                                    % x and y are vectors
% Convert vectors to matrices
[X, Y] = meshgrid(x, y);
% Flip upside down the Y matrix because of the tif file...
Y = flipud(Y);

%generate ice_elevation:

h220 = dlmread('thickness_380y.txt');
topo = dlmread('domain.dat');
ice_elevation = topo + h220;              
contourf(ice_elevation,60);
colorbar

%interpolate points for glacier surface profile:

flowline_elevation = interp2(X,Y,ice_elevation,rx,ry)

%interpolate points for bed surface:

bed_elevation = interp2(X,Y,topo,rx,ry)

%modify the flowline_elevation array from a row to a column: 

flowline_elevation = flowline_elevation';



%%%%%%%%%%%%%% plot surface ice velocity (m/yr) only for a specific glacier by cutting polygon and calculate max, minimum and mean surface velocity statistics for given glacier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This code extracts the ice thickness output data to include only the data
%within a polygon imported as a shapefile, which was previously drawn in
%ArcGIS to highlight a specific area: currently the El Loro glacier catchment
%All that needs to be set is the location and name of the Domain (tiffile),
%the location and name of the polygon (shapefile), and the location and
%name of the surfspeed files (files) and its correct folder (study
%directory). The code then computes surface ice velocity for each surfspeed files
%present in the folder and also produces a PNG image for each of them.

shapefile = 'I:\tanc\tanc_field\Model_runs\Cordillera Blanca\Cojup_valley\Holocene_simulations\Shapefiles\glacier_polygon\Cojup_glacier_polygon.shp';
tiffile = 'I:\tanc\tanc_field\domains\ALOS world 3d DEM\Cordillera_blanca\Cojup_glacier\WGS84\30m\no_modern_glacier\Lake_bump_trimmed\lake_removed_manually\larger_dem_for_holocene_simulations\LIA_moraine_removed\Moraine_removed_DEM.tif';

polygon = shaperead(shapefile);
% When read, polygon.X, polygon.Y have trailing NaN, we need to get rid of
% it
rx = polygon.X(1: end-1); % take all array elements except the last one (NaN)
% same for Y:
ry = polygon.Y(1: end-1);

% Read tif file
geo = geotiffinfo(tiffile);
% Get the pixel coordinates (convert pixels to X,Y)
[x,y] = pixcenters(geo);                                                               % x and y are vectors
% Convert vectors to matrices
[X, Y] = meshgrid(x, y);
% Flip upside down the Y matrix because of the tif file...
Y = flipud(Y);

% Generate the mask: confront the polygon with the tiffile to get what is
% inside the polygon and what is outside
mask = inpolygon(X,Y, rx, ry);                                                         % mask is a logical (1,0) matrix
mask = double(mask);                                                                   % convert logical to double so as to multiply it with the elevation

% Read all surfspeed files
study_directory = workenvironment;
files = dir(strcat(study_directory, 'surfSpeed_*y.txt'));
volume = [];

for k=1:length(files)
    % Read the surfSpeed file
    h0 = dlmread(strcat(study_directory, files(k).name));

    %masked h0:
    res = times(h0, mask);                                                               % multiply each cell by its equivalent cell in other matrix
    res(res>400) = 0;                                                                    % get rid of anomalies (crasy high values): this value in bracket (m/s) might change from one glacier to the next
    res(res == 0) = NaN;                                                                 % get rid of 0 values in order to only visualize the glacier extent as well as velocity
    pcolor(res); shading flat; colorbar
    title(['surface velocity (m/yr)']);
    colormap jet;
    highest_volicity(k) = max(res(:));                                                   % calculates maximum ice velocity
    lowest_velocity(k) = min(res(:));                                                    % calculates minimum ice velocity
    Mean_velocity(k) = mean(res, 'all','omitnan');                                       % calculates mean surface ice velocity across the glacier, thus for the entire polygon-cut matrix ('all'), and by omitting the NaN values ('omitnan').
    filename = strrep(files(k).name, 'surfSpeed_', 'Cojup_glacier_surface_velocity_');   % Replace in the filename 'thickness_' by 'Rio_comissario_only_' (nothing)
    filename = strsplit(filename,'.txt');                                                % remove .txt extension from filename
    title(filename, 'Interpreter', 'none');                                              % Interpreter: none means that the _ will not be interpreted as a subscript
    filenamepng = strcat(filename{1}, '.png');
    saveas(gcf, filenamepng)                                                             % export current figure (gca) into a .png file rather than a .fig file
    clf % clear figure ?
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%