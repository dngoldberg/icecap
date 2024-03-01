# icecap
a simple, isothermal shallow-ice matlab mddel for fast modelling of mountain glaciers
Shallow Ice Cap Model
Instructions and basic MATLAB guidance
D Goldberg, 19 May 2017

RECENTLY COPIED FROM WORD, NEED TO EDIT

The ice cap flow model is written in the MATLAB programming language. MATLAB is in many ways similar to Python, and also in many ways easier to use – it simply makes use of slightly different conventions and syntax so it might take time to get used to.

RUNNING THE MODEL
The ice model is run by calling the function icecap which is in the file icecap.m. At the prompt, type 
 icecap (1)
The reason for typing (1) means that you are specifying runmode=1 – an ice model run (alternative described below). This will begin running the model, and making plots, as you have seen before. It will also create data files, with a frequency (in years) given by freq_files in params.m. The files will have the names
thickness_Ny.txt
fluxE_Ny.txt
fluxW_Ny.txt
fluxS_Ny.txt
fluxN_Ny.txt
Where N is the time in the model in years. They are data (text) files that can be read into MATLAB arrays similarly to how you have done in Python. The contents of the Thickness file should be clear. The others attempt to capture velocity (or rather flux) information. The model grid is divided into cells (small boxes), the number and size of which depend on (a) the size of your input topography file and (b) the grid size parameters Lx and Ly. The size of each cells is (x, y) where x is given by Lx divided by the number of east-west pixels in your topography file. fluxE is a map of the mass flux per unit width (velocity times thickness, units m2/year) through the east through all of these cells. Note that flux is positive if it is out of the cell. Thus if ice is moving eastward, you will see positive values for fluxE but negative for fluxW. The fluxes  for a given cell are shown here:
 


params.m contains a number of run parameters, which I have tried to explain through comments (green text). If any are non intuitive let me know. Remember that each time you change params.m, you must save the file before running the model, or the changes will not be picked up.
To end a model prematurely at any point, press Ctrl-C.
Try running the model with the default paramters. Change the length of the run to 50 years in params.m. Now try reading a thickness file and plotting the result. At the command prompt type the following commands:
>> icecap(1)
time stepping model
>> h50 = dlmread('thickness_40y.txt');
>> contourf(h50);
>> colorbar
The first command runs the model. The second defines a 2D array h50 and reads into it the contents of thickness_40y.txt. The next command makes a filled contour plot, and the final command adds a colorbar. If you prefer a different color scheme, type “colormap jet” at the next prompt. What do you see? Try reading one of the flux files and do the same thing. 

SURFACE MASS BALANCE
The ice model implements 3 different parameterisations of surface mass balance: DEFAULT, KO (Kaser and Osmaston), and DD (Degree Day). Which one you use is chosen by the parameter mbal_type (1=Default, 2=KO, 3=DD). 
DEFAULT: this was previously the only option. Here mass balance at elevation z is given by 
b = –  T, T>0
b = P,        T<0
T = (zELA – z) × /1000 + TELA,
TELA = P / 
where T is temperature, P is precipitation,  is the melt factor, and  is lapse rate is deg C/km. Note that TELA is found from the parameters, and Tsea_level is not given. The parameters P, zELA, and are specified in the DEFAULT section of params.m.  applies to all mass balance choices.
KO: Here a theoretical profile is constructed without precipitation. zELA is calculated from sea level temperature and lapse rate, and mass balance has an elevation-mass balance slope in the ablation zone, and a different slope in the accumulation zone. These two slopes are specified in the KO section of params.m, while sea level temperature is specified with Tsl in params.m.
DD: This specification is very similar to DEFAULT. 
b = P –  T, T>0
b = P,           T<0
T = (zELA) × /1000 + TSL,
Note here, there is precipitation where T>0. Although the formulation uses the same parameters, the parameters relevant to DD appear in a separate section. That is, melt_factor and DD_melt_factor can be set independently; the former is only relevant when mbal_type=1, and the latter when mbal_type=3.
You can view a mass balance – elevation profile without running the simulation. This can be done by calling
>> icecap(0,zmin,zmax)
Where zmin and zmax need to be given as bounds for the elevation plot (you will get an error otherwise).

BASAL SLIDING
The ice model may be specified as frozen bed (zero sliding velocity) or with a uniform sliding parameter of the form
us = Asd p
where us is sliding velocity (which is additive with deformational velocity), As is a constant parameter (set by A_weert in params.m; if absent from the file it is set to zero by default), d is driving stress, and p is an exponent specified by A_weert in params.m (generally set to 3). As = 10-15 and p=3 gives a sliding velocity of 100 m/a when driving stress is 1 bar (100 kPa).

TIMESTEPPING AND INSTABILITY
The glacier thickness evolution is a diffusion equation, meaning it is subject to numerical instability if the time step parameter (dt in params.m) is too large. This can result in nonphysical features (see below the “stripes” in flux which correspond to nonphysical oscillations in thickness), or negative or infinite thickness values; at this point the simulation has diverged from a representation of glacier physics. The ice-sheet equations are solved implicitly which means the velocity which drives thickness change and the new thickness are solved simultaneously, extending the valid time step range somewhat, but the case shown below can still happen.
There will not be simple relationship that will predict how small the time step should be; it is related to thickness and velocity, and will be different depending on how much velocity is due to sliding. Empirically I find a good timestep is 1/5 of the time it takes for the highest velocity to move across a pixel, but a bit more investigation is needed (including the investigation of numerical parameters which are beyond the scope of this manual). In the simulation shown below, pixel size is 60m, maximum speed in about 300m/a, and the time step is 0.05a, which is too large. The instability does not occur when time step is 0.025a.



PLOTTING YOUR RESULTS
If you have completed a model run and would like to re-make the plots that are produced during the run, you can do this with
>> icecap(2,time)
Where time should be a value (in years) for which you know files have been saved (e.g. thickness_10y.txt). You will get an error if the file does not exist, or if you do not provide the time argument.

ARRAY INDEXING
Arrays are less confusing than in Python (I think) because indices start with 1 instead of zero. It is also important to note that for a 2D array, the row is always specified first. 
If I want to know the value of row 47, column 50 of the thickness array I read in earlier, I would type

>> h50(47,50)
79.9050

If I wanted to see the entire 47th row I would type
>> h50(47,:)
But this would lead to a lot of numbers of the screen and I could not use them later. Rather I can define a new 1D arrays
>> hrow47 = h50(47,:);
(end the commands with a semicolon if you want to decrease output.) Then if I want to take a look at these numbers visually I can type
>> plot (hrow47); 
What happens when you do this?
Another task you might need to do is to know how large an array is. Type 
>> size(h50) 
(with no semicolon). You will see two numbers: the number of rows and the number of columns, in that order. What is the size of hrow47?

SOME BASIC PROGRAMMING ELEMENTS
For-loops are very similar to Python. In the program forLoop.m I have written a simple code that runs a loop, which defines an array and plots it. In fact it does the same thing twice, but using different approaches – the second is a way to “concatenate” arrays. Run it, look through it, try to make sense of it! One element of this module is the command “hold on”. This allows multiple lines on one plot. Without this, each “plot” command erases the old plot. “hold off” takes this away. In python, you had to explicitly refresh a plot otherwise it would draw over what was already there.
Note: when running modules, any variable you define will still be around after you run the module (see the “workspace” tab in the console). To clear all the variables in the workspace, type “clear all”. Python does not allow this.

Now, once you are comfortable with this, here are some exercises  you can do.
1.  Try reading the topography file given to you into an array, and (a) create a filled contour plot and (b) create a line plot of the 50th row. This is simple enough to do on the command line.
2. Create a module which will read all of the thickness files created by the module in sequence, and create a 1D array which holds the total ice volume at each time level, and plot the array. You will start with a for loop – how will the steps be defined? To read a file in the for loop, you will need to construct an appropriate file name. This can be done by concatenating a string as forLoop concatenated an array, eg
filename = [‘thickness_’ num2str(30) ‘y.txt’];
the above converts “30” to a string which is placed in the larger string; you can then call 
dlmread(filename) 
and use this to read the appropriate data into an array. Of course, you will need to define your for loop to do this for every appropriate time level. Next, to turn a thickness array H into a volume, you can call
volume = 1000 * 1000* sum(sum(H))
the 1000 * 1000 factor is the area of a cell (at least with the current settings). “sum” is a function that finds the sum of an array. The inner call to “sum” finds the sum of each column of H, creating a 1D array. The outer call to “sum” finds the sum of this array. 
If you are successful, try with longer model runs (you will need to modify the code module you created). How long (how many years) before the volume is steady?




