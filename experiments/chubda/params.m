% MODEL RUN PARAMETERS

% number of years of the simulation
nYears = 200;    
start_year = 1900;

% time step (years) -- time in between updates of ice cap velocity and
% surface
dt = 0.25;   

% how often you would like thickness and flux files written (years)
freq_files = 10;

% how often you would like to plot the ice thickness etc
freq_plot = 10;

% the softness parameter in the ice flow law, Pa^(-3) y^(-1)
Aglen = 1e-16;    

% the exponent in the ice flow law
nglen = 3;

% the weertman sliding constant and exponent; fits into equation
% u_sliding = A_s taud^p, has units m / (Pa^n-year)

p_weert = 3;
A_weert = 1e-15; % gives sliding velocity of 100m/a where tau = 1 bar


% DAMAGE PARAMETERS
use_damage_param = true;
damage_factor = .1; % (m/a)^-1: fraction of damage per meter per year of negative mass balance
damage_max = 0.6;    % maximum damage fraction

% ice density 
rho = 917;
% ice density x gravity, N/m^3
rhog = rho * 9.8;

output_speed_to_file = true; % or false!


% file which will be read to define topography -- given in meters
topofile = 'domain.dat';
% topofile = 'testTopo';

% file which will be read to define ocean location -- 1 where there is
% land, 0 where there is ocean
oceanmaskfile = 'mask.dat';
% oceanmaskfile = 'testMask';

% width of pixel in x-direction (m)
dx = 50.0302;  % resolution of the DEM provided by OGGM
% Lx=300*300;

% width of pixel in y-direction (m)
dy = 50.0302;  % resolution of the DEM provided by OGGM
% Ly = 300*300;

% optional: x coordinate of bottom left pixel
x_global = 2.717590438725613e+05;
% optional: y coordinate of bottom left pixel
y_global = 3.099184185798137e+06;


% determine if initial thickness to be read from file
%thickness_path = 'thickness_2090y.txt';

% MASS BALANCE PARAMETERS

% choice of mass balance parameterisation.
% 1 = DEFAULT: ELA specified. See documentation for details.
% 2 = KO: Kaser and Osmaston (see Reay dissertation, Fig 2.1)
% 3 = DD: Degree day representation, very similar to Default. see doc.

mbal_type = 3; 

% atmospheric lapse rate (deg C per km) -- needed for all models
lapse_rate = 5.5;
% sea level temp-- needed for models 2 & 3 -- deg C
Tsl = 35.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for DEFAULT
%%%%%%%%%%%%%%%%%%%%%%%%

% ELA elevation (m)
elev_ELA    = 5000;

clim_from_file = true;
clim_file = 'climatedata.csv';


% melt-temp factor (m/a / deg C) melt rate per deg C above zero
melt_factor = 0.31 ;

% constant accumulation rate (m/a)
accum = 1.15 ;

%%%%%%%%%%%%%%%%%%%%
% Parameters for KO%
%%%%%%%%%%%%%%%%%%%%

% accumulation zone MB profile -- in meters elev per (kg/m^2/yr)
KO_slope_acc = 834;
% ablation zone MB profile -- in meters elev per (kg/m^2/yr)
KO_slope_abl = 50;

%%%%%%%%%%%%%%%%%%%%
% Parameters for DD%
%%%%%%%%%%%%%%%%%%%%

% melt factor: given in meters ice per yr per degree C. Only causes melt
% when T>0
DD_melt_factor = .29;
% precipitation: given in meters ice per yr
DD_precip = 1.2;



plot_smb=true;
