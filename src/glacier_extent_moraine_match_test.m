%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load bed_elevation table (if already in workspace) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If workspace already contains bed_elevation_table.mat

%load the table: 

load('bed_elevation_table.mat')

%then run the surface gradient to get flowline elevation code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shallap glacier surface gradient code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Open shapefile:

shapefile = 'C:\Users\gg1tpl\Desktop\Tancrede_leger\PHD\Patagonia\Ice_Sheet_Modelling\Ice_cap_modelling\tanc_field\Model_runs\Cordillera_Blanca\Cojup_valley\Holocene_simulations\Shapefiles\centreline\Cojup_centreline.shp';
tiffile = 'C:\Users\gg1tpl\Desktop\Tancrede_leger\Peru-CASCADA-project\Icecap_modelling_CASCADA\AW3D30_DEM\Cojup\Final_DEM_after_modern_glacier_and_lake_removal\DEM_for_early_holocene_simulations\LIA_moraine_removal\point_to_DEM_IDW\Moraine_removed_DEM.tif';

flowline = shaperead(shapefile);

%get rid of NaN

rx = flowline.X(1: end-1);
ry = flowline.Y(1: end-1);

geo = geotiffinfo(tiffile);                                                     % Read tif file
[x,y] = pixcenters(geo); % x and y are vectors                                  % Get the pixel coordinates (convert pixels to X,Y)
[X, Y] = meshgrid(x, y);                                                        % Convert vectors to matrices                                                                         
Y = flipud(Y);                                                                  % Flip upside down the Y matrix because of the tif file...

                                                                                %generate ice_elevatio
h260 = dlmread('thickness_380y.txt');
topo = dlmread('domain.dat');
ice_elevation = topo + h260;              
contourf(ice_elevation,60);
colorbar

flowline_elevation = interp2(X,Y,ice_elevation,rx,ry)                           %interpolate points for glacier surface profile:
bed_elevation = interp2(X,Y,topo,rx,ry)                                         %interpolate points for bed surface:

                                                                                %Add data to table

flowline_elevation = flowline_elevation';                                       %modify the flowline_elevation array from a row to a column: 
flowline_elevation_table = array2table(flowline_elevation)                      %convert the flowline_elevation double array to a table: 
Table_bed_elevation = [Table_bed_elevation flowline_elevation_table];           %add the new flowline table (1column) to the table_bed_elevation table. 
glacier_surface_gradient_table = Table_bed_elevation;                           %Copy to new variable with new and different name.
save('glacier_surface_gradient_table.mat', 'glacier_surface_gradient_table')    %Save the new version of the Table with the new column: new name


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Plot glacier elevation linegraph  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('glacier_surface_gradient_table.mat')                                      %Load the right table
x1 = glacier_surface_gradient_table{:,2};                                       %Set the X1 variable as 2nd table column: i.e distance to first point in (m). 
x2 = glacier_surface_gradient_table{:,4};                                       %Set the x2 variable as 4th table column: i.e. the target moraine distance to headwall 
y1 = glacier_surface_gradient_table{:,6};                                       %Set the Y1 variable as 6th table column: i.e. flowline elevation
y2 = glacier_surface_gradient_table{:,3};                                       %Set the Y2 variable as 3rd table column: i.e. bed elevation
y3 = glacier_surface_gradient_table{:,5};                                       %Set the Y3 variable as 5th table column: i.e. moraine elevation (display a line for moraine location)
plot(x1,y1,x1,y2,x2,y3);                                                        %plot a linegraph and decide on settings: i.e. title etc
title('glacier surface elevation');
xlabel('distance to headwall (m)');
ylabel('Elevation (m)');
xlim([5800,6400]);
ylim([4250,4650]);