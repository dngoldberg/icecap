function M = SMB (s_elev)

global elev_ELA melt_factor lapse_rate accum time

    T_ela = accum/melt_factor;
    theta = (elev_ELA-s_elev) * lapse_rate/1000 + T_ela;
    M = accum-melt_factor * theta;
    M(theta<0) = accum;

    
    
    %M = mbal_slope * (s_elev - elev_ELA);
    %M (s_elev>elev_thresh) = mbal_slope * (elev_thresh - elev_ELA);
    
return
    