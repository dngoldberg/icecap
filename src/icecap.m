function icecap (runMode,arg1,arg2,arg3)

% runMode = 1: run model
% runMode = 0: 

if (runMode==1)
    disp('time stepping model');
elseif (runMode==0)
    if (nargin==3)
     min_elev_disp = arg1;
     max_elev_disp = arg2;
     time_in=-9999999;
    elseif (nargin==4)
     min_elev_disp = arg1;
     max_elev_disp = arg2;
     time_in = arg3;
    else
     error('provide zmin and zmax for display of SMB profile in m: icecap(0,zmin,zmax)');
    end
elseif (runMode==2)
    if (nargin==2)
     time_for_plot = arg1;
    else
     error('provide time for plot in years: icecap(0,time)');
    end
end
   

% load parameters

%ensure parameters will be available to other functions
global topofile oceanmaskfile mbal_type lapse_rate Tsl elev_ELA ...
 melt_factor accum KO_slope_acc KO_slope_abl DD_melt_factor ...
 DD_precip DD_melt_temp rho time clim_from_file slt_series acc_series clim_yr

params;

% check parameters

if (dt < 1)
    dt = 1 / round(1/dt);
else
    dt = round(dt);
end

disp(['timestep set to' num2str(dt) ' years']);

freq_plot = round(freq_plot/dt) * dt;
if (freq_plot == 0) 
    error('freq_plot must be a multiple of dt');
end

freq_files = round(freq_files/dt) * dt;
if (freq_files == 0) 
    error('freq_files must be a multiple of dt');
end


C = 2*Aglen/(nglen+2)*(rhog)^nglen;
Csurf = 2*Aglen/(nglen+1)*(rhog)^nglen;
if(~exist('p_weert','var') | ~exist('A_weert','var'));
    C_slid = 0;
else
    if (~exist('p_weert','var'))
        p_weert = 3;
    end
    C_slid = (rhog)^p_weert * A_weert;
end

%initialisation of topo and grid

[topo mask nx ny] = readSpatialData;
Lx = nx * dx;
Ly = ny * dy;
if(~exist('x_global','var'));
    x_global = 0;
end
if(~exist('y_global','var'));
    y_global = 0;
end
X = x_global+linspace(dx/2,Lx-dx/2,nx);
Y = y_global+linspace(dy/2,Lx-dy/2,ny);
[X Y] = meshgrid(X,Y);

dlmwrite('Xgrid.txt',X);
dlmwrite('Ygrid.txt',Y);


% if climate file, initialisation of arrays

if (clim_from_file);

 MM = csvread(clim_file);
 slt_series = MM(:,2);
 acc_series = MM(:,3);
 clim_yr = MM(:,1);
 [xx ii] = sort(clim_yr);
 slt_series = slt_series(ii);
 acc_series = acc_series(ii);

end


% ADVANCED -- set up arrays for numerical solver


tic

if (runMode==1)

disp('preparing model arrays');    
    
H = zeros(size(topo));


if (exist('thickness_path','var'));
    if (length(thickness_path)>0);
        ss=dir(thickness_path);
        if (length(ss)>0);
            disp(['reading initial thickness ' thickness_path]);
            Hinit = dlmread(thickness_path);
            if (size(Hinit)==size(H))
                H = Hinit;
            else
                error('initial thickness wrong size');
            end
        else
            disp(['initial thickness file ' thickness_path ' not found']);
        end
    end
end


Hexp = zeros(ny+2,nx+2);
topoExp = zeros(ny+2,nx+2);
topoExp(2:end-1,2:end-1) = topo;

flow2West = zeros(size(H));
flow2East = zeros(size(H));
flow2North = zeros(size(H));
flow2South = zeros(size(H));


for i=1:ny
    
    for j=1:nx
        if (j>1)
            flow2West(i,j) = mask(i,j)==1 & mask(i,j-1)==1;
        end
        if (j<nx)
            flow2East(i,j) = mask(i,j)==1 & mask(i,j+1)==1;
        end
        if (i>1)
            flow2South(i,j) = mask(i,j)==1 & mask(i-1,j)==1;
        end
        if (i<ny)
            flow2North(i,j) = mask(i,j)==1 & mask(i+1,j)==1;
        end
    end
end
flow2West = logical(flow2West);
flow2East = logical(flow2East);
flow2North = logical(flow2North);
flow2South = logical(flow2South);


end

% find continentality factor at every point on grid

%landFac = land_factor(X,Y,mask); 
landFac = ones(size(X));

tt = toc
disp(['model arrays prepared, ' num2str(tt) ' seconds. beginning stepping']);


%timestepping of ice thickness

if (runMode==1)
    
time = 0;
timestep = 0;

n_inner_loop_max = 10;
inner_loop_tol = 1e-6;


if (exist('start_year'));
 time = start_year;
end

while (timestep < nYears*round(1/dt))
    
% ADVANCED -- set up arrays for numerical solver
  switch(mbal_type);
        case 1
            M = SMB(H + topo);
        case 2
            M = SMB_KO(H + topo);
        case 3
            M = SMB_DD(H + topo);
  end

  H_inner = H;
  
  Aglen_array = Aglen * ones(size(H));
  
  if (exist('use_damage_param') & exist('damage_factor') & exist('damage_max'))
      if(use_damage_param)
          damage_frac = zeros(size(H));
          damage_frac(M<0) = -damage_factor * M(M<0);
          damage_frac(damage_frac>damage_max) = damage_max;
          Aglen_array(M<0) = Aglen_array(M<0) ./ (1-damage_frac(M<0));
      end
  end
  
  Carray = 2*Aglen_array/(nglen+2)*(rhog)^nglen;
  Csurfarray = 2*Aglen_array/(nglen+1)*(rhog)^nglen;
     
  for k_inner = 1:n_inner_loop_max;
      
    H_last_iter = H_inner;
    
    Hmid = (H_inner+H)/2;   
    
    Hexp(2:end-1,2:end-1) = Hmid;   
    
    Sexp = topoExp + Hexp;
    
    D = zeros(size(Hexp));
	Dsurf = zeros(size(Hexp));
    Sx = (Sexp(2:end-1,3:end)-Sexp(2:end-1,1:end-2))/2/dx;
    Sy = (Sexp(3:end,2:end-1)-Sexp(1:end-2,2:end-1))/2/dy;
    D(2:end-1,2:end-1) = ...
          Carray .* (Sx.^2 + Sy.^2+1e-8) .^ ((nglen-1)/2) .* Hmid.^(nglen+2);
    D(2:end-1,2:end-1) = D(2:end-1,2:end-1) + ...
          C_slid * (Sx.^2 + Sy.^2+1e-8) .^ ((nglen-1)/2) .* Hmid.^(nglen+1);

    Dsurf(2:end-1,2:end-1) = ...			
	      Csurfarray .* (Sx.^2 + Sy.^2+1e-8) .^ ((nglen-1)/2) .* Hmid.^(nglen+1) + ...
		  C_slid * (Sx.^2 + Sy.^2+1e-8) .^ ((nglen-1)/2) .* Hmid.^(nglen);
		  
    Dx = .5 * (D(2:end-1,1:end-1)+D(2:end-1,2:end));
    Dy = .5 * (D(1:end-1,2:end-1)+D(2:end,2:end-1));
	Dx_surf = .5 * (Dsurf(2:end-1,1:end-1)+Dsurf(2:end-1,2:end));
    Dy_surf = .5 * (Dsurf(1:end-1,2:end-1)+Dsurf(2:end,2:end-1));

    DxWest = Dx(:,1:end-1);
    DxEast = Dx(:,2:end);
    DyNorth = Dy(2:end,:);
    DySouth = Dy(1:end-1,:);
	
	DxWest_surf = Dx_surf(:,1:end-1);
    DxEast_surf = Dx_surf(:,2:end);
    DyNorth_surf = Dy_surf(2:end,:);
    DySouth_surf = Dy_surf(1:end-1,:);

    A00 = ones(size(H));
    A00(mask==1) = A00(mask==1)  ...
        + dt/dx^2 * DxWest(mask==1) ...
        + dt/dx^2 * DxEast(mask==1) ...
        + dt/dy^2 * DyNorth(mask==1) ...
        + dt/dy^2 * DySouth(mask==1);

    Ap1 = zeros(size(H));
    Ap1 (flow2North) = Ap1(flow2North) - dt/dy^2 * DyNorth(flow2North);

    Am1 = zeros(size(H));
    Am1 (flow2South) = Am1(flow2South) - dt/dy^2 * DySouth(flow2South);

    ApN = zeros(size(H));
    ApN (flow2East) = ApN(flow2East) - dt/dx^2 * DxEast(flow2East);

    AmN = zeros(size(H));
    AmN (flow2West) = AmN(flow2West) - dt/dx^2 * DxWest(flow2West);

    A00 = reshape(A00,nx*ny,1);
    Ap1 = reshape(Ap1,nx*ny,1);
    Am1 = reshape(Am1,nx*ny,1);
    ApN = reshape(ApN,nx*ny,1);
    AmN = reshape(AmN,nx*ny,1);

    AA = spdiags([ApN Ap1 A00 Am1 AmN],[-ny -1 0 1 ny],nx*ny,nx*ny)';


    rhs = reshape(H,nx*ny,1);

    rhs = rhs  + ...
         reshape( DxEast .* (topoExp(2:end-1,3:end)-topoExp(2:end-1,2:end-1)) * dt/dx^2 ...
                - DxWest .* (topoExp(2:end-1,2:end-1)-topoExp(2:end-1,1:end-2)) * dt/dx^2 ...
                + DyNorth .* (topoExp(3:end,2:end-1)-topoExp(2:end-1,2:end-1)) * dt/dy^2 ...
                - DySouth .* (topoExp(2:end-1,2:end-1)-topoExp(1:end-2,2:end-1)) * dt/dy^2, ...
        nx*ny,1);

    H_inner = reshape(AA \ rhs,ny,nx);
    neg_vol = sum(H(H<0)) * dx * dy;
    if (neg_vol>0)
        disp(['Volume lost: ' num2str(-neg_vol) ' cubic m']);
    end
    H_inner(H_inner<0)=0;
    
    err_iter = norm(H_inner-H_last_iter,2)/nx/ny;
    %disp(num2str(err_iter))

    if (err_iter < inner_loop_tol)
        %disp(['convergence of inner loop: ' num2str(k_inner) ' iterations']);

        break
           
    end
    
  end
  if(k_inner==n_inner_loop_max)
      %disp(['did not converge: ' num2str(err_iter) ' error']);      
  end
  
H = H_inner;

H = H + M * dt;

H(H<0) = 0;


fluxNorth = zeros(size(H));
fluxSouth = zeros(size(H));
fluxWest  = zeros(size(H));
fluxEast  = zeros(size(H));

velNorth = zeros(size(H));
velSouth = zeros(size(H));
velWest  = zeros(size(H));
velEast  = zeros(size(H));


S = H + topo;

filecheck_time = mod(timestep,round(freq_files/dt));
plot_time = mod(timestep,round(freq_plot/dt));

if (filecheck_time == 0 | plot_time == 0)

    Hexp(2:end-1,2:end-1) = H;    
    Sexp = topoExp + Hexp;
	
    fluxNorth = -((Sexp(3:end,2:end-1)-Sexp(2:end-1,2:end-1))/dy) .* DyNorth;
    fluxSouth = ((Sexp(2:end-1,2:end-1)-Sexp(1:end-2,2:end-1))/dy) .* DySouth;

    fluxEast = -((Sexp(2:end-1,3:end)-Sexp(2:end-1,2:end-1))/dx) .* DxEast;
    fluxWest = ((Sexp(2:end-1,2:end-1)-Sexp(2:end-1,1:end-2))/dx) .* DxWest;
	
	velNorth = -((Sexp(3:end,2:end-1)-Sexp(2:end-1,2:end-1))/dy) .* DyNorth_surf;
    velSouth = -((Sexp(2:end-1,2:end-1)-Sexp(1:end-2,2:end-1))/dy) .* DySouth_surf;

    velEast = -((Sexp(2:end-1,3:end)-Sexp(2:end-1,2:end-1))/dx) .* DxEast_surf;
    velWest = -((Sexp(2:end-1,2:end-1)-Sexp(2:end-1,1:end-2))/dx) .* DxWest_surf;
	    
end

vel_u = .5 * velEast + .5 * velWest;
vel_v = .5 * velNorth + .5 * velSouth;

output_speed = sqrt(vel_u.^2 + vel_v.^2);

if (filecheck_time == 0)

    filename_end = ['_' num2str(round(time)) 'y.txt'];
    dlmwrite(['thickness' filename_end],H);
    dlmwrite(['fluxN' filename_end],fluxNorth);
    dlmwrite(['fluxS' filename_end],fluxSouth);
    dlmwrite(['fluxE' filename_end],fluxEast);
	dlmwrite(['fluxW' filename_end],fluxWest);
	if (output_speed_to_file)
	 dlmwrite(['surfSpeed' filename_end],output_speed);
	end
end




% END OF THICKNESS UPDATE


% PLOT THICKNESS

if (plot_time==0)
    
    figure(1)
    make_plot (time, H, topo, output_speed, [], M);
    
    
end;


timestep = timestep + 1;
time = start_year + timestep*dt;
disp(['time: ' num2str(time) ' yrs']);
end

elseif(runMode==0)
    
    close (figure(1))
    
    time = time_in;
    Z = linspace(min_elev_disp,max_elev_disp,1000);
    switch(mbal_type)
        case 1
            M = SMB(Z);
        case 2
            M = SMB_KO(Z);
        case 3
            M = SMB_DD(Z);
    end
    
    plot(M,Z); grid on; ylabel('elevation');
    set(gca,'xlim',[min(M) max(M)*1.2]);
    xlabel(['SMB']);
    if (time_in>0);
        title (['SMB in year ' num2str(time)]);
    else
        title('SMB')
    end
else
    
    time = time_for_plot;
    filename_end = ['_' num2str(round(time)) 'y.txt'];
    s=dir(['thickness' filename_end]);
    if (length(s)~=1)
     error('output file not found: please note the "freq_files" parameter value and choose a suitable time');
    end
    H = dlmread(['thickness' filename_end]);
    fluxEast = dlmread(['fluxE' filename_end]);
    fluxSouth = dlmread(['fluxS' filename_end]);
    readSpatialData;
    Z = H + topo;
    switch(mbal_type)
        case 1
            M = SMB(Z);
        case 2
            M = SMB_KO(Z);
        case 3
            M = SMB_DD(Z);
    end
    figure()
    make_plot (time, H, topo, fluxEast, fluxSouth, M);

end

return
