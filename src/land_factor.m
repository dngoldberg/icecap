function landFac = land_factor(X,Y,mask); 

global contRadius

 landFac = zeros(size(mask));
 [ny nx] = size(mask);
 for i=1:ny;
     
     for j=1:nx;
         
         x = X(i,j);
         y = Y(i,j);
         supp = ( sqrt(((X-x).^2 + (Y-y).^2)) < contRadius );
         areaLand = sum(mask(supp)==1);
         area = sum(sum(supp));
         landFac(i,j) = areaLand/area;
     end
 end
return

