
% Script to make movie of velocity contours and level set evolution.

function movies(start,numframes)

close all
warning off all

%load 'centers.dat'

% Open an avi file

%aviobjx = avifile('x_cont_movie.avi') 
aviobjy = avifile('ycontours_pluslevelset2.avi'); 
%aviobjp = avifile('intermediate.avi') 
%aviobjpress = avifile('press.avi') 

cirx = [-.05:.001:.05];
ciry = sqrt(.0025-cirx.^2);

filenum = start;

for j = start:numframes
    
        filenum = filenum +1  
    
    if filenum < 10

        filestring = ['file00', int2str(filenum),'.dat'];
        
    elseif filenum < 100
        
        filestring = ['file0', int2str(filenum),'.dat'];

    else
        
        filestring = ['file', int2str(filenum),'.dat'];
        
    end
    
    fid = fopen(filestring,'r');
    
    xnodes = str2num(fgetl(fid));
    ynodes = str2num(fgetl(fid));
    
    for i = 1:100
        INT(i) = str2double(fgetl(fid)); 
	end
    
    totalpts = xnodes*ynodes;

    data = fscanf(fid,'%g',[15,totalpts]);
    
    X = data(1,:)';
    Y = data(2,:)';
    PHI = data(3,:)';
    S_PHI = data(4,:)';
    W_PHI = data(5,:)';
    H = data(6,:)';
    S_H = data(7,:)';
    W_H = data(8,:)';
    U = data(10,:);
    V = data(11,:);
    PHIX = data(12,:)';
    PHIY = data(13,:)';
    PHIXY = data(14,:)';
        
    
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
        PHIX_for_plot(i,:) = PHIX(xnodes*(i-1)+1:xnodes*i);
        PHIY_for_plot(i,:) = PHIY(xnodes*(i-1)+1:xnodes*i);
        PHIXY_for_plot(i,:) = PHIXY(xnodes*(i-1)+1:xnodes*i);
    end 
    
%     hx =    figure(1)
%     contour(X_for_plot,Y_for_plot,U_for_plot,[2 2],'-m')
%     hold on
%     contour(X_for_plot,Y_for_plot,U_for_plot,[16/9 16/9],'--m')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[14/9 14/9],':m')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[12/9 12/9],'-r')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[10/9 10/9],'--r')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[8/9 8/9],':r')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[6/9 6/9],'-y')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[4/9 4/9],'--y')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[2/9 2/9],':y')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[0 0],'-g')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-1 -1],'-k')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-8/9 -8/9],'--k')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-7/9 -7/9],':k')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-6/9 -6/9],'-b')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-5/9 -5/9],'--b')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-4/9 -4/9],':b')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-3/9 -3/9],'-c')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-2/9 -2/9],'--c')
%     contour(X_for_plot,Y_for_plot,U_for_plot,[-1/9 -1/9],':c')    
    
      hy =    figure(1);
      set(hy, 'Position', [800 0 1000 1000])
%     contour(X_for_plot,Y_for_plot,V_for_plot,[9/18 9/18],'-m')
%     hold on
%     contour(X_for_plot,Y_for_plot,V_for_plot,[8/18 8/18],'--m')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[7/18 7/18],':m')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[6/18 6/18],'-r')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[5/18 5/18],'--r')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[4/18 4/18],':r')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[3/18 3/18],'-y')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[2/18 2/18],'--y')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[1/18 1/18],':y')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[0 0],'-g')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-1 -1],'-k')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-8/18 -8/18],'--k')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-7/18 -7/18],':k')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-6/18 -6/18],'-b')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-5/18 -5/18],'--b')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-4/18 -4/18],':b')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-3/18 -3/18],'-c')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-2/18 -2/18],'--c')
%     contour(X_for_plot,Y_for_plot,V_for_plot,[-1/18 -1/18],':c')
    
     hpress =    figure(1);
      set(hpress, 'Position', [800 0 1000 1000])
 %     contour(X_for_plot,Y_for_plot,U_for_plot,20)
      %quiver((X_for_plot(:,2:end)+X_for_plot(:,1:end-1))./2,Y_for_plot(:,1:end-1),U_for_plot(:,1:end-1),zeros(size(U_for_plot(:,1:end-1))))
      hold on
      %quiver(X_for_plot(1:end-1,:),(Y_for_plot(2:end,:)+Y_for_plot(1:end-1,:))./2,zeros(size(V_for_plot(1:end-1,:))),V_for_plot(1:end-1,:))
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[45 45],'-m')
 %    hold on
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[40 40],'--m')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[35 35],':m')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[30 30],'-r')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[25 25],'--r')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[20 20],':r')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[15 15],'-y')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[10 10],'--y')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[5 5],':y')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[0 0],'-g')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-5 -5],'-k')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-10 -10],'--k')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-15 -15],':k')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-20 -20],'-b')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-25 -25],'--b')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-30 -30],':b')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-35 -35],'-c')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-40 -40],'--c')
 %    contour(X_for_plot,Y_for_plot,P_for_plot,[-45 -45],':c')


 %    hp =    figure(3)
    %hold on
    %surf(X_for_plot,Y_for_plot,PHI_for_plot) 
    hold on
     contour(X_for_plot,Y_for_plot,PHI_for_plot,[0 0],'k')
    %hold on
    % contour(X_for_plot,Y_for_plot,S_PHI_for_plot,[0 0],'k')
    %hold on
    contour(X_for_plot,Y_for_plot,W_PHI_for_plot,[0 0],'k')
    hold on
     %surf(X_for_plot,Y_for_plot,1.0-H_for_plot)
     %quiver(X_for_plot,Y_for_plot,U_for_plot,V_for_plot)
     %quiver(X_for_plot,Y_for_plot,PHIX_for_plot,PHIY_for_plot)

    if INT(2) == 1 
        cgx = INT(50);
        cgy = INT(51);
        th = INT(52);
        rot_mat = [cos(th) -sin(th); ...
                   sin(th) cos(th)];
        internal = rot_mat*[0 .05]';
        plot([cgx cgx+internal(1)], [cgy cgy+internal(2)], '-k', ...
             cgx + cirx, cgy + ciry, 'r', cgx + cirx, cgy - ciry, 'r')

    elseif INT(2) == 2 | INT(2) == 3
	plot(INT(50), INT(51),'.', ...
	     INT(52), INT(53),'o')
    
    end
    %axis equal
    axis([-0.1 9.2 -0.1 0.85])
    %axis([-0.1 1.1 -0.1 1.1])
 %%   hy =    figure(1)
      hpress = figure(1);
%     diameter = 0.4;
%     w = diameter;
%     h = diameter;
%     [a b] = find(S_PHI_for_plot == max(max(S_PHI_for_plot)))
%     x = centers(filenum+1,1) + 0.3;
%     y = centers(filenum+1,2) + 0.4;
% %    x = 0.4;
% %    y = 0.4;
%     rectangle('Position',[x y w h],'Curvature',[1,1],'FaceColor','k') 
       
    %aviobjx = addframe(aviobjx,hx)
    aviobjy = addframe(aviobjy,hpress);
    %aviobjp = addframe(aviobjp,hp)
    %aviobjpress = addframe(aviobjpress,hpress)
    fclose(fid);
    clf
end

close all
%aviobjx = close(aviobjx)
aviobjy = close(aviobjy)
%aviobjp = close(aviobjp)
%aviobjp = close(aviobjpress

end
