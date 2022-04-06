function M = SMB_DD (s_elev)

global lapse_rate DD_melt_factor DD_precip Tsl time ...
 clim_from_file slt_series acc_series clim_yr

    prec = DD_precip;
	temp = Tsl;
 
    if (clim_from_file);
	 if (time >= min(clim_yr) & time < max(clim_yr));
	  yr0_ind = max(find(clim_yr<=time));
	  yr1_ind = min(find(clim_yr>time));
	  temp = (time-clim_yr(yr0_ind)) * slt_series(yr1_ind) + ...
	        (clim_yr(yr1_ind)-time) * slt_series(yr0_ind);
	  prec = (time-clim_yr(yr0_ind)) * acc_series(yr1_ind) + ...
	        (clim_yr(yr1_ind)-time) * acc_series(yr0_ind);
	 end
	end
	  

    T_surf = temp*ones(size(s_elev)) - lapse_rate * s_elev/1000;
    M = zeros(size(s_elev));

    M(T_surf>0) = -DD_melt_factor * T_surf(T_surf>0);
    M = ...
        M + prec;    
    
return
