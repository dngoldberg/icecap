function [topo mask nx ny] = readSpatialData

global topofile oceanmaskfile

    topo = dlmread(topofile);
    mask = dlmread(oceanmaskfile);
    topo(mask==0) = 0;
    [ny nx] = size(topo);
    
return

