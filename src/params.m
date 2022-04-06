% MODEL RUN PARAMETERS

% number of years of the simulation
nYears = 1000;     

% time step (years) -- time in between updates of ice cap velocity and
% surface
dt = 0.1 * 1/2;   

% how often you would like thickness and flux files written (years)
freq_files = 50;

% how often you would like to plot the ice thickness etc
freq_plot = 5;

% decide if you would like to output surface speed
output_speed_to_file = true;

% the softness parameter in the ice flow law, Pa^(-3) y^(-1)
Aglen = 1e-16;    

% the exponent in the ice flow law
nglen = 3;

% the weertman sliding constant and exponent; fits into equation
% u_sliding = A_s taud^p, has units m / (Pa^n-year)

p_weert = 3;
A_weert = 1e-15; % gives sliding velocity of 100m/a where tau = 1 bar

% ice density 
rho = 917;
% ice density x gravity, N/m^3
rhog = rho * 9.8;




% file which will be read to define topography -- given in meters
topofile = 'domain.dat';
% topofile = 'testTopo';

% file which will be read to define ocean location -- 1 where there is
% land, 0 where there is ocean
oceanmaskfile = 'mask.dat';
% oceanmaskfile = 'testMask';

% width of pixel in x-direction (m)
dx = 60;
% Lx=300*300;

% width of pixel in y-direction (m)
dy = 60;
% Ly = 300*300;

% optional: x coordinate of bottom left pixel
x_global = 296496.975591494;
% optional: y coordinate of bottom left pixel
y_global = 5148034.7188424;


% MASS BALANCE PARAMETERS

% choice of mass balance parameterisation.
% 1 = DEFAULT: ELA specified. See documentation for details.
% 2 = KO: Kaser and Osmaston (see Reay dissertation, Fig 2.1)
% 3 = DD: Degree day representation, very similar to Default. see doc.

mbal_type = 3; 

% atmospheric lapse rate (deg C per km) -- needed for all models
lapse_rate = 6.5;
% sea level temp-- needed for models 2 & 3 -- deg C
Tsl = 20;


% the following, if "true", will override any setting 
%  of temperature and accum via parameters during
%  the period defined within the file
% IMPORTANT -- only implemented with Degree day representation
clim_from_file = true;
% file should give annual temperature and accumulation
% year in 1st column, temperature in 2nd column
% annual accumulation in 3rd column, no headers
% should also set "start_year" below
clim_file = 'Book1.csv';
start_year = 2020;


%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for DEFAULT
%%%%%%%%%%%%%%%%%%%%%%%%

% ELA elevation (m)
elev_ELA    = 1650;

% melt-temp factor (m/a / deg C) melt rate per deg C above zero
melt_factor = 0.3 ;

% constant accumulation rate (m/a)
accum = 8.0 ;

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
DD_melt_factor = .5;
% precipitation: given in meters ice per yr
DD_precip = 5.27;




