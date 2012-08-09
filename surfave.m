% Script to make movie of velocity contours and level set evolution.

clear
close all

load 'centers.dat'
% ave = []
% 
% xnodes = 46;
% ynodes = 46;
% 
% totalpts = xnodes*ynodes;
% 
% filenum = -1
% 
% for j = 1:201
%     
%     filenum = filenum +1;
%     
%     if filenum < 10
% 
%         filestring = ['file00', int2str(filenum),'.dat'];
%         
%     elseif filenum < 100
%         
%         filestring = ['file0', int2str(filenum),'.dat'];
%         
%     else
%         
%         filestring = ['file', int2str(filenum),'.dat'];
%     end
%     
%     fid = fopen(filestring,'r');
%     
%     line1 = fgetl(fid)
%     line2 = fgetl(fid)
% 
%     data = fscanf(fid,'%g',[9,totalpts]);
%     
%     X = data(1,:)';
%     Y = data(2,:)';
%     PHI = data(3,:)';
% %    S_PHI = data(4,:)';
% %    P = data(7,:)';
% %    U = data(8,:)';
% %    V = data(9,:)';
%     
%     for i = 1:ynodes
%         X_for_plot(i,:) = X(xnodes*(i-1)+1:xnodes*i);
%         Y_for_plot(i,:) = Y(xnodes*(i-1)+1:xnodes*i);
%         PHI_for_plot(i,:) = PHI(xnodes*(i-1)+1:xnodes*i);
%  %       S_PHI_for_plot(i,:) = S_PHI(xnodes*(i-1)+1:xnodes*i);
%   %      U_for_plot(i,:) = U(xnodes*(i-1)+1:xnodes*i);
%   %      V_for_plot(i,:) = V(xnodes*(i-1)+1:xnodes*i);
%   %      P_for_plot(i,:) = P(xnodes*(i-1)+1:xnodes*i);
%     end 
%     
%     [c,h] = contour(X_for_plot,Y_for_plot,PHI_for_plot,[0 0],'k');
%     
%     temp = size(c,2)
%     
%     surf_vals = c(2,2:temp);
%     pos_vals = c(1,2:temp);
%     clipped_val = [];
%     
%     for m = 1:temp-1
%         if (pos_vals(m) < 0.3) || (pos_vals(m) > 0.7)
%             clipped_val = [clipped_val surf_vals(m)]
%         end
%     end
%     
%     ave(j) = mean(clipped_val);
%     
%     fclose(fid)
% 
%     close all
% end

x = [0.005:0.005:4];
y = centers(:,2);
y = y + 0.6;
%y = y - ave';
plot(x,y,x,0.4848,'k',x,0.4835,'g',x,0.4520,'r');
%hold on
%plot(x,0.4848,'k',x,0.4835,'g',x,0.4520,'r')
%plot(x,0.4835,'g')
%plot(x,0.4520,'r')
legend('numerical solutions','analytical buoyancy solution','adjusted for fluid leakage','adjusted for loss of fluid volume')

    