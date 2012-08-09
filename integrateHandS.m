% Script for MATLAB postprocessing data that comes out of the boiling code
% Specifically for integrating the volume of fluid.  
% MPX
% Jessica Sanders
% 06/16/2008

clear
close all

big_area = 0
little_area = 0

xnodes = 31;
ynodes = 31;

delx = 1/(xnodes-1);
dely = 1/(ynodes-1);

totalpts = xnodes*ynodes;

    fid = fopen('file500.dat','r');

    line1 = fgetl(fid)
    line2 = fgetl(fid)

    data = fscanf(fid,'%g',[9,totalpts]);

    X = data(1,:)';
    Y = data(2,:)';
    PHI = data(3,:)';
    S_PHI = data(4,:)';
    H = data(5,:)';
    S_H = data(6,:)';
    P = data(7,:)';
    U = data(8,:)';
    V = data(9,:)';
    
    for i = 1:ynodes
        X_for_plot(i,:) = X(xnodes*(i-1)+1:xnodes*i);
        Y_for_plot(i,:) = Y(xnodes*(i-1)+1:xnodes*i);
        PHI_for_plot(i,:) = PHI(xnodes*(i-1)+1:xnodes*i);
        S_PHI_for_plot(i,:) = S_PHI(xnodes*(i-1)+1:xnodes*i);      
        H_for_plot(i,:) = H(xnodes*(i-1)+1:xnodes*i);
        S_H_for_plot(i,:) = S_H(xnodes*(i-1)+1:xnodes*i);        
        U_for_plot(i,:) = U(xnodes*(i-1)+1:xnodes*i);
        V_for_plot(i,:) = V(xnodes*(i-1)+1:xnodes*i);
        P_for_plot(i,:) = P(xnodes*(i-1)+1:xnodes*i);
    end 
    
    for i = 1:(xnodes-1)
        for j = 1:(ynodes-1)
    
        % Get the value of H at the center point.
    
        aveH = 0.25*(H_for_plot(i,j) + H_for_plot(i+1,j) + H_for_plot(i+1,j+1) + H_for_plot(i,j+1));
        aveSH = 0.25*(S_H_for_plot(i,j) + S_H_for_plot(i+1,j) + S_H_for_plot(i+1,j+1) + S_H_for_plot(i,j+1));
    
        %multipliy by the area
    
        little_area = aveH*delx*dely - aveH*aveSH*delx*dely;
    
        big_area = big_area + little_area;
    
        end
    end
    
    big_area
    
