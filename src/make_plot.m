function make_plot (time, H, topo, speed, dummy, mBal);

    params;
    figure(1);

    [ny nx] = size(topo);
    Lx = nx * dx;
    Ly = ny * dy;
    X = linspace(dx/2,Lx-dx/2,nx);
    Y = linspace(dy/2,Lx-dy/2,ny);
    [X Y] = meshgrid(X/1000,Y/1000);
    
    if (exist('use_damage_param') & exist('damage_factor') & exist('damage_max'))
      if(use_damage_param)
          damage_frac = zeros(size(H));
          damage_frac(mBal<0) = -damage_factor * mBal(mBal<0);
          damage_frac(damage_frac>damage_max) = damage_max;
          damage_frac(H==0) = nan;
          
      end
    end
    
    Hplot = H;
    Hplot(H==0)=nan;
    subplot(2,3,1); 
    pcolor(X,Y,Hplot+topo); shading flat;
    xlabel('X (km)');
    ylabel('Y (km)');
    title(['elevation, ' num2str(time) ' years']);
    colorbar

    colormap jet;
    hold on
    contour(X,Y,topo,5,'k','linewidth',2);
    hold off
    
    subplot(2,3,2); 

if ~(isempty(dummy));    
    fluxEast = speed;
    fluxSouth = dummy;
    
    xvel = zeros(size(fluxEast));
    xvel(:,1:end-1) = 2*fluxEast(:,1:end-1)./(H(:,1:end-1)+H(:,2:end));
    yvel = zeros(size(fluxSouth));
    yvel(2:end,:) = 2*fluxSouth(2:end,:)./(H(1:end-1,:)+H(2:end,:));
%     
    xvel(isnan(Hplot))=nan;
    yvel(isnan(Hplot))=nan;
    
    maxc = min( max(max(abs(xvel(~isnan(Hplot)))),max(abs(yvel(~isnan(Hplot))))), 200);
    
    pcolor(X,Y,xvel); shading flat; colorbar
    title(['xvel (m/a), ' num2str(time) ' years']);
    caxis([-maxc maxc]); 
    
    
    subplot(2,3,3);
    pcolor(X,Y,-yvel); shading flat; colorbar
    title(['yvel (m/a), ' num2str(time) ' years']);
    caxis([-maxc maxc]); 
    
else
    speed(isnan(Hplot))=nan;
    pcolor(X,Y,speed); shading flat; colorbar
    title(['speed (m/a), ' num2str(time) ' years']);
    maxc = min( max(speed(~isnan(Hplot))),200);
    caxis([0 maxc]); 
    
    
end

    
    
    subplot(2,3,4);
    pcolor(X,Y,mBal); shading flat; colorbar; 
    title(['surf mass balance, ' num2str(time) ' years']);
    


    subplot(2,3,4);
    pcolor(X,Y,Hplot); shading flat; colorbar
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
    if (exist('plot_smb'));
     if(plot_smb)
      pcolor(X,Y,M0); shading flat;
      title('Surface Mass Balance (m/yr)')
      colorbar; 
     else
      pcolor(X,Y,topo); shading flat; colorbar
      title(['topography'])
      hold on;
      contour(M0,[0 0],'k','linewidth',2);
      hold off;
     end
    else
     pcolor(X,Y,topo); shading flat; colorbar
     title(['topography'])
     hold on;
     contour(M0,[0 0],'k','linewidth',2);
     hold off;
    end
    
    
    
    if (time>start_year)
        last_file_time = floor((time-dt)/freq_files)*freq_files;
        h_last = dlmread(['thickness_' num2str(last_file_time) 'y.txt']);
        subplot(2,3,6);
        last_file_time
        time-last_file_time
        pcolor((H-h_last)/(time-last_file_time)); shading flat; colorbar
        title(['rate of change (m/a)']);
    end
    
    if (exist('use_damage_param') & exist('damage_factor') & exist('damage_max'))
      if(use_damage_param)
          figure(2);
          pcolor(X,Y,damage_frac);
          xlabel('X (km)');
          ylabel('Y (km)');
          title('damage fraction')
          shading flat
          colorbar
      end
    end
    
    drawnow
