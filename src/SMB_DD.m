function M = SMB_DD (s_elev)

global lapse_rate DD_melt_factor DD_precip DD_melt_temp Tsl

    T_surf = Tsl*ones(size(s_elev)) - lapse_rate * s_elev/1000;
    M = zeros(size(s_elev));

    M(T_surf>0) = -DD_melt_factor * T_surf(T_surf>0);
    M = ...
        M + DD_precip;
    
    
    
return