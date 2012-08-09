
close all

filestring = 'bench.dat';
    
fid = fopen(filestring,'r');

data = fscanf(fid,'%g',[5,200000]);

time = data(1,:)';
ind = data(2,:)';
inc = data(3,:)';
pos = data(4,:)';
vol = data(5,:)';

fclose(fid);

%texp = [0.43 0.62 0.8 0.97 1.14 1.29 1.45 1.62 1.76 1.93 2.07 2.24 2.4 2.54 2.71 2.87 3.04 3.21 3.29 3.32]';
%zexp = [1.11 1.22 1.44 1.67 1.89 2.11 2.33 2.56 2.78 3 3.22 3.44 3.67 3.89 4.11 4.33 4.56 4.78 4.89 5.0]';

%plot(time/0.07633,pos/0.05715,texp,zexp,'o')

plot(time,vol/vol(2))