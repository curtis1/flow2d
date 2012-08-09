% Script for MATLAB postprocessing data that comes out of the boiling code
% MPX
% Jessica Sanders
% 06/16/2008

clear
close all

xnodes = 76;
ynodes = 16;

totalpts = xnodes*ynodes;

%filenum = 29;

%for j = 1:15

%    filenum = filenum +1;
%    filestring = ['file0', int2str(filenum),'.dat'];
    
%    fid = fopen(filestring,'r');
    
    fid = fopen('file190.dat','r');

    line1 = fgetl(fid)
    line2 = fgetl(fid)

    data = fscanf(fid,'%g',[12,totalpts]);

    X = data(1,:)';
    Y = data(2,:)';
    PHI = data(3,:)';
    S_PHI = data(4,:)';
    W_PHI = data(5,:)';
    H = data(6,:)';
    S_H = data(7,:)';
    W_H = data(8,:)';
    P = data(9,:)';
    U = data(10,:)';
    V = data(11,:)';
    
    for i = 1:ynodes
        X_for_plot(i,:) = X(xnodes*(i-1)+1:xnodes*i);
        Y_for_plot(i,:) = Y(xnodes*(i-1)+1:xnodes*i);
        PHI_for_plot(i,:) = PHI(xnodes*(i-1)+1:xnodes*i);
        S_PHI_for_plot(i,:) = S_PHI(xnodes*(i-1)+1:xnodes*i); 
        W_PHI_for_plot(i,:) = W_PHI(xnodes*(i-1)+1:xnodes*i); 
        H_for_plot(i,:) = H(xnodes*(i-1)+1:xnodes*i);
        S_H_for_plot(i,:) = S_H(xnodes*(i-1)+1:xnodes*i);
        W_H_for_plot(i,:) = W_H(xnodes*(i-1)+1:xnodes*i); 
        U_for_plot(i,:) = U(xnodes*(i-1)+1:xnodes*i);
        V_for_plot(i,:) = V(xnodes*(i-1)+1:xnodes*i);
        P_for_plot(i,:) = P(xnodes*(i-1)+1:xnodes*i);
    end 

%    subplot(3,5,j)
    figure(1)
    contour(X_for_plot,Y_for_plot,PHI_for_plot,[0 0],'k','LineWidth',1)
    %figure(1)
    %surf(X_for_plot,Y_for_plot,PHI_for_plot)
    hold on
    contour(X_for_plot,Y_for_plot,S_PHI_for_plot,[0 0],'k','LineWidth',1)
    contour(X_for_plot,Y_for_plot,W_PHI_for_plot,[0 0],'k','LineWidth',1)
    %figure(2)
    %surf(X_for_plot,Y_for_plot,S_PHI_for_plot)
    axis equal
    axis([0 5 0 1])
%end






