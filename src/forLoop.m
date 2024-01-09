% one way to create the array y

for i=5:5:40;
    
    disp(i)
    j = round(i/5);  % always a good idea to round when numbers will be used for indices
    y(j) = i^2;
    
end


% another way to create the same array    


y2 = [];
for i=5:5:40;
    
    disp(i)
    y2 = [y2 i^2];
    
end

plot(y,'b');
hold on;
plot(y2,'r+');

