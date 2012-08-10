
% Script to make movie of velocity contours and level set evolution.

function aframe(start)

%close all
warning off all

cirx = [-.05:.001:.05];
ciry = sqrt(.0025-cirx.^2);

filenum = start;
 
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
P = data(9,:)';
U = data(10,:);
V = data(11,:);
PHIX = data(12,:)';
PHIY = data(13,:)';
PHIXY = data(14,:)';
FSIG = data(15,:)';


for i = 1:ynodes
    X_for_plot(i,:) = X(xnodes*(i-1)+1:xnodes*i);
    Y_for_plot(i,:) = Y(xnodes*(i-1)+1:xnodes*i);
    PHI_for_plot(i,:) = PHI(xnodes*(i-1)+1:xnodes*i);
    S_PHI_for_plot(i,:) = S_PHI(xnodes*(i-1)+1:xnodes*i); 
    W_PHI_for_plot(i,:) = W_PHI(xnodes*(i-1)+1:xnodes*i); 
    H_for_plot(i,:) = H(xnodes*(i-1)+1:xnodes*i);
    S_H_for_plot(i,:) = S_H(xnodes*(i-1)+1:xnodes*i); 
    W_H_for_plot(i,:) = W_H(xnodes*(i-1)+1:xnodes*i); 
    P_for_plot(i,:) = P(xnodes*(i-1)+1:xnodes*i); 
    U_for_plot(i,:) = U(xnodes*(i-1)+1:xnodes*i);
    V_for_plot(i,:) = V(xnodes*(i-1)+1:xnodes*i);
    PHIX_for_plot(i,:) = PHIX(xnodes*(i-1)+1:xnodes*i);
    PHIY_for_plot(i,:) = PHIY(xnodes*(i-1)+1:xnodes*i);
    PHIXY_for_plot(i,:) = PHIXY(xnodes*(i-1)+1:xnodes*i);
    FSIG_for_plot(i,:) = FSIG(xnodes*(i-1)+1:xnodes*i);
end 

  hy =    figure(1);
%  set(hy, 'Position', [800 0 1000 1000])

hpress =    figure(1);
%set(hpress, 'Position', [800 0 1000 1000])
%     contour(X_for_plot,Y_for_plot,U_for_plot,20)
%quiver((X_for_plot(:,2:end)+X_for_plot(:,1:end-1))./2,Y_for_plot(:,1:end-1),U_for_plot(:,1:end-1),zeros(size(U_for_plot(:,1:end-1))))
hold on
%quiver(X_for_plot(1:end-1,:),(Y_for_plot(2:end,:)+Y_for_plot(1:end-1,:))./2,zeros(size(V_for_plot(1:end-1,:))),V_for_plot(1:end-1,:))

%    hp =    figure(3)
%hold on
%surf(X_for_plot,Y_for_plot,H_for_plot) 
hold on
contour(X_for_plot,Y_for_plot,PHI_for_plot,100)
hold on
contour(X_for_plot,Y_for_plot,PHI_for_plot,[0 0],'k')
% contour(X_for_plot,Y_for_plot,S_PHI_for_plot,[0 0],'k')
%hold on
 contour(X_for_plot,Y_for_plot,W_PHI_for_plot,[0 0],'k')
hold on
%surf(X_for_plot,Y_for_plot,sqrt(PHIX_for_plot.^2+PHIY_for_plot.^2))
%surf(X_for_plot,Y_for_plot,FSIG_for_plot)
%surf(X_for_plot(4:end-4,4:end-4),Y_for_plot(4:end-4,4:end-4),FSIG_for_plot(4:end-4,4:end-4))
%quiver(X_for_plot,Y_for_plot,U_for_plot,V_for_plot)
quiver(X_for_plot,Y_for_plot,PHIX_for_plot,PHIY_for_plot)

%if INT(2) == 1 
%    cgx = INT(50);
%    cgy = INT(51);
%    th = INT(52);
%    rot_mat = [cos(th) -sin(th); ...
%               sin(th) cos(th)];
%    internal = rot_mat*[0 .05]';
%    plot([cgx cgx+internal(1)], [cgy cgy+internal(2)], '-k', ...
%         cgx + cirx, cgy + ciry, 'r', cgx + cirx, cgy - ciry, 'r')

%elseif INT(2) == 2 | INT(2) == 3
%plot(INT(50), INT(51),'.', ...
%     INT(52), INT(53),'o')
%
%end

axis equal
axis([-0.000 9.12 -0.000 0.76])
%axis([-0.000 0.00625 -0.000 0.0125])
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
%aviobjp = addframe(aviobjp,hp)
%aviobjpress = addframe(aviobjpress,hpress)
fclose(fid);

%aviobjp = close(aviobjpress

end
