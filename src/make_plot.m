function make_plot (time, H, topo, fluxEast, fluxSouth, mBal);

    params;

    [ny nx] = size(topo);
    Lx = nx * dx;
    Ly = ny * dy;
    X = linspace(dx/2,Lx-dx/2,nx);
    Y = linspace(dy/2,Lx-dy/2,ny);
    [X Y] = meshgrid(X,Y);
    
    Hplot = H;
    Hplot(H==0)=nan;
    subplot(2,3,1); 
    pcolor(X,Y,Hplot+topo); shading flat;
    xlabel('X (m)');
    ylabel('Y (m)');
    title(['elevation, ' num2str(time) ' years']);
    colorbar

    colormap jet;
    hold on
    contour(topo,'k','linewidth',2);
    hold off
    subplot(2,3,2); 

    xvel = zeros(size(fluxEast));
    xvel(:,1:end-1) = 2*fluxEast(:,1:end-1)./(H(:,1:end-1)+H(:,2:end));
    yvel = zeros(size(fluxSouth));
    yvel(2:end,:) = 2*fluxSouth(2:end,:)./(H(1:end-1,:)+H(2:end,:));
    
    maxc = max(max(max(abs(xvel),abs(yvel))));
    
    pcolor(xvel); shading flat; colorbar
    title(['x-velocity (m/a), ' num2str(time) ' years']);
    clim("auto"); 
    
    subplot(2,3,3);
    
    pcolor(-yvel); shading flat; colorbar
    title(['y-velocity (m/a), ' num2str(time) ' years']);
    clim("auto"); 
    
    %subplot(2,3,4);
    %pcolor(mBal); shading flat; colorbar; 
    %title(['surf mass balance, ' num2str(time) ' years']);


    subplot(2,3,4);
    pcolor(Hplot); shading flat; colorbar
    title(['thickness, ' num2str(time) ' years']);
    colormap jet
   
    switch(mbal_type);
        case 1
            M0 = SMB(topo);
        case 2
            M0 = SMB_KO(topo);
        case 3
            M0 = SMB_DD(topo);
    end
    
    subplot(2,3,5); 
    pcolor(topo); shading flat; colorbar
    title(['topography'])
    hold on;
    contour(M0,[0 0],'k','linewidth',2);
    hold off;
    
    
    
    if (time>0)
        last_file_time = floor((time-dt)/freq_files)*freq_files;
        h_last = dlmread(['thickness_' num2str(last_file_time) 'y.txt']);
        subplot(2,3,6);
        last_file_time
        time-last_file_time
        pcolor((H-h_last)/(time-last_file_time)); shading flat; colorbar
        title(['rate of change (m/a)']);
    end
    
    drawnow
