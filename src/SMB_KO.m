function M = SMB_KO (s_elev)

global lapse_rate KO_slope_acc KO_slope_abl Tsl rho

    T_surf = Tsl*ones(size(s_elev)) - lapse_rate * s_elev/1000;
    z_ela = Tsl/lapse_rate*1000;
    M = zeros(size(s_elev));
    M(s_elev<z_ela) = ...
        (s_elev(s_elev<z_ela)-z_ela)/KO_slope_abl;
    M(s_elev>z_ela) = ...
        (s_elev(s_elev>z_ela)-z_ela)/KO_slope_acc;
    
return