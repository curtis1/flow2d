% Curtis Lee
% 2/7/11

function lsgint(xn,yn,deltaT,endtime,plotstyle)  % Input dimensions (mesh size) & duration


% Domain Variables

hx = 1/xn;                      % Horizontal spacing
hy = 1/yn;                      % Vertical spacing
x = 0:hx:1;                     % x-nodes
y = 0:hy:1;                     % y-nodes
X = ones(yn+1,1)*x;             % x-node matrix
Y = y'*ones(1,xn+1);            % y-node matrix
u = zeros(size(X));             % Horizontal velocity
v = zeros(size(X));             % Vertical velocity
T = 8;                          % Velocity-field time constant
aviobjy = avifile('ls.avi'); 

% For refined contour imaging

if (plotstyle == 1 || plotstyle == 3 || plotstyle == 5)
    refine = 3;
    xr = xn.*refine;
    yr = yn.*refine;
    XR = ones(yr+1,1)*[0:1/xr:1];
    YR = [0:1/yr:1]'*ones(1,xr+1);
    phifine = zeros(size(XR));
end
% Velocity in x & y

ufun = @(x,y,time) - 2.*cos(pi.*time./T).*(sin(pi.*x).^2).*sin(pi.*y).*cos(pi.*y);
vfun = @(x,y,time)   2.*cos(pi.*time./T).*(sin(pi.*y).^2).*sin(pi.*x).*cos(pi.*x);

u = ufun(X,Y,0);
v = vfun(X,Y,0);

% Deformation Matrix

vmat1 = @(x,y,time) -    pi.*cos(pi.*time./T).*sin(2.*pi.*x).*sin(2.*pi.*y);
vmat2 = @(x,y,time)   2.*pi.*cos(pi.*time./T).*(sin(pi.*y).^2).*cos(2.*pi.*x);
vmat3 = @(x,y,time) - 2.*pi.*cos(pi.*time./T).*(sin(pi.*x).^2).*cos(2.*pi.*y);
vmat4 = @(x,y,time)      pi.*cos(pi.*time./T).*sin(2.*pi.*x).*sin(2.*pi.*y);

% Update Variables

%deltaT = 0.4.*max(max([u v])).*min(hx, hy)./2

X1 = [0 0];                     % Space Update 1
X2 = [0 0];                     % Space Update 2
X0 = [0 0];                     % Space Update Final

I = [1 0; 0 1];                 % Identity

gX1 = [0 0; 0 0];               % Grad Update 1
gX2 = [0 0; 0 0];               % Grad Update 2
gX0 = [0 0; 0 0];               % Grad Update Final

phin1 = zeros(size(X));         % NEW phi value

phigradn1 = [0; 0];             % NEW gradient at a single point
phigxn1 = zeros(size(X));       % NEW x-dir gradient
phigyn1 = zeros(size(X));       % NEW y-dir gradient

% Initial & Gradient Level-set parameters

r0 = 0.15;  % 0.15
x0 = 0.5;   % 0.5
y0 = 0.75;  % 0.75

phi0 = @(x,y) exp(-(x-x0).^2-(y-y0).^2) - exp(-r0.^2);

phigx0 = @(x,y) - 2.*(x - x0).*exp(-(x-x0).^2-(y-y0).^2);
phigy0 = @(x,y) - 2.*(y - y0).*exp(-(x-x0).^2-(y-y0).^2);
phigxy0 = @(x,y) 4.*(x - x0).*(y - y0).*exp(-(x-x0).^2-(y-y0).^2);


phi = phi0(X,Y);
phigx = phigx0(X,Y);
phigy = phigy0(X,Y);
phigxy = phigxy0(X,Y);

if (plotstyle == 1 || plotstyle == 3)
    phifine0 = phi0(XR,YR);
    phifine = phifine0;
end
if (plotstyle == 5)
    phigxi = phigx0(XR,YR);
    phigyi = phigy0(XR,YR);
end
    
% Plot Initial Configuration
    
figure(1)
set(1,'OuterPosition',[0 500 1200 1000])

if plotstyle == 1        % Refined Cubic Interpolation (Better Contour)
    contour(XR,YR,phifine,[0 0],'w','LineWidth',2)
    surface(XR,YR,phifine)
    view([-.7,-.7,1]);
elseif plotstyle == 2    % View using only calculated points of phi (Poor contour)
    contour(X,Y,phi,[0 0],'w','LineWidth',2)
    surface(X,Y,phi)
    view([-.7,-.7,1]);
elseif plotstyle == 3    % 2D with refined contour
    %quiver(x,y,phigx, phigy)
    %hold on
    contourf(XR,YR,phifine0,[0 0],'g','LineWidth',2)
    hold on
    contour(XR,YR,phifine,[0 0],'k','LineWidth',2)
elseif plotstyle == 4
    quiver(x,y,phigx, phigy)
    hold on
    contour(x,y,phigxy)
    hold on
    contour(X,Y,phi,[0 0],'r','LineWidth',2)
elseif plotstyle == 5
    quiver(XR,YR,phigxi, phigyi)
    hold on
    contour(x,y,phigxy)
    hold on
    contour(X,Y,phi,[0 0],'r','LineWidth',2)
else                % 2D with poor contour
    quiver(x,y,phigx, phigy)
    hold on
    contour(X,Y,phi,[0 0],'w','LineWidth',2)
end


hold on
axis equal
axis([0 1 0 1])
pic = figure(1);
aviobjy = addframe(aviobjy,pic);
figure(1)
%clf(1)

% BEGIN TIME LOOP
tic
for time = deltaT:deltaT:endtime
toc
tic
    % Find Cross-Derivatives (Cell-based approach)  <--- THIS DIDN'T WORK!
    
 %   for j = 1:yn
 %       for i = 1:xn
 %           
 %           bottom = (phigy(j,i+1) - phigy(j,i))./hx;
 %           top = (phigy(j+1,i+1) - phigy(j+1,i))./hx;
 %           left = (phigx(j+1,i) - phigx(j,i))./hy;
 %           right = (phigx(j+1,i+1) - phigx(j,i+1))./hy;
 %           
 %           phigxy(j,i) = 0.75.*(left + bottom) - 0.25.*(top + right);
 %           
 %       end
 %   end
 %   
 %   
 %   for j = 1:yn
 %       
 %       i = xn+1;
 %       
 %       bottom = (phigy(j,i) - phigy(j,i-1))./hx;
 %       top = (phigy(j+1,i) - phigy(j+1,i-1))./hx;
 %       left = (phigx(j+1,i-1) - phigx(j,i-1))./hy;
 %       right = (phigx(j+1,i) - phigx(j,i))./hy;
 %           
 %       phigxy(j,i) = 0.75.*(right + bottom) - 0.25.*(top + left);
 %           
 %   end
 %   
 %   for i = 1:xn
 %       
 %       j = yn+1;
 %       
 %       bottom = (phigy(j-1,i+1) - phigy(j-1,i))./hx;
 %       top = (phigy(j,i+1) - phigy(j,i))./hx;
 %       left = (phigx(j,i) - phigx(j-1,i))./hy;
 %       right = (phigx(j,i+1) - phigx(j-1,i+1))./hy;
 %           
 %       phigxy(j,i) = 0.75.*(left + top) - 0.25.*(bottom + right);
 %           
 %   end
 %   
 %   i = xn+1;
 %   j = yn+1;
 %   
 %   bottom = (phigy(j-1,i) - phigy(j-1,i-1))./hx;
 %   top = (phigy(j,i) - phigy(j,i-1))./hx;
 %   left = (phigx(j,i-1) - phigx(j-1,i-1))./hy;
 %   right = (phigx(j,i) - phigx(j-1,i))./hy;
 %           
 %   phigxy(j,i) = 0.75.*(right + top) - 0.25.*(bottom + left);
    
    % Find Cross Derivatives (Central Difference Method)  <--- WORKS!
    
    for j = 1:yn+1
        for i = 2:xn
            
            phigxy(j,i) = (phigy(j,i+1) - phigy(j,i-1))./(2.*hx);
            
        end
    end
    
    for j = 2:yn
        
        phigxy(j,1) = (phigx(j+1,1) - phigx(j-1,1))./(2.*hy);
        phigxy(j,yn+1) = (phigx(j+1,yn+1) - phigx(j-1,yn+1))./(2.*hy);
                
    end
     
    phigxy(1,1) = (phigxy(1,2) + phigxy(2,1))./2;           % These can be more accurate.
    phigxy(xn+1,1) = (phigxy(xn,1) + phigxy(xn+1,2))./2;
    phigxy(1,yn+1) = (phigxy(1,yn) + phigxy(2,yn+1))./2;
    phigxy(xn+1,yn+1) = (phigxy(xn+1,yn) + phigxy(xn,yn+1))./2;
     
    % Calculate velocity at grid points
    
    %u = ufun(X,Y,time);
    %v = vfun(X,Y,time);
    








    % Superconsistent Shu-Osher RK3 Scheme
    
    for j = 1:yn+1
        for i = 1:xn+1
            
            X1 = [X(j,i) Y(j,i)] - deltaT.*[Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,1) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,4)];
            
            gX1 = I - deltaT.*[Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,2) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,5); ...
                               Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,3) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,6)];
                           
            X2 = [X(j,i) Y(j,i)] - 0.25.*deltaT.*[Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,1) + Hi(hx,hy,X,Y,X1(1),X1(2),time,1) ...
                                                  Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,4) + Hi(hx,hy,X,Y,X1(1),X1(2),time,4)];
                                              
            gX2 = I - 0.25.*deltaT.*([Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,2) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,5); ...
                                      Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,3) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,6)] ...
                               + gX1*[Hi(hx,hy,X,Y,X1(1),X1(2),time,2) Hi(hx,hy,X,Y,X1(1),X1(2),time,5); ...
                                      Hi(hx,hy,X,Y,X1(1),X1(2),time,3) Hi(hx,hy,X,Y,X1(1),X1(2),time,6)]);
                                                            
            X0 = [X(j,i) Y(j,i)] - (1./6).*deltaT.*[Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,1) + Hi(hx,hy,X,Y,X1(1),X1(2),time,1) ...
                                                                    + 4.*Hi(hx,hy,X,Y,X2(1),X2(2),time+0.5.*deltaT,1) ...
                                                    Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,4) + Hi(hx,hy,X,Y,X1(1),X1(2),time,4) ...
                                                                    + 4.*Hi(hx,hy,X,Y,X2(1),X2(2),time+0.5.*deltaT,4)];
                                                                
            gX0 = I - (1./6).*deltaT.*([Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,2) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,5); ...
                                      Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,3) Hi(hx,hy,X,Y,X(j,i),Y(j,i),time+deltaT,6)] ...
                               + gX1*[Hi(hx,hy,X,Y,X1(1),X1(2),time,2) Hi(hx,hy,X,Y,X1(1),X1(2),time,5); ...
                                      Hi(hx,hy,X,Y,X1(1),X1(2),time,3) Hi(hx,hy,X,Y,X1(1),X1(2),time,6)] ...
                            + 4.*gX2*[Hi(hx,hy,X,Y,X2(1),X2(2),time+deltaT/2,2) Hi(hx,hy,X,Y,X2(1),X2(2),time+deltaT/2,5); ...
                                      Hi(hx,hy,X,Y,X2(1),X2(2),time+deltaT/2,3) Hi(hx,hy,X,Y,X2(1),X2(2),time+deltaT/2,6)]);                      
                                  
            phin1(j,i) = H(hx,hy,X,Y,phi,phigx,phigy,phigxy,j,i,X0(1),X0(2),1);
            
            phigradn1 = gX0 * [H(hx,hy,X,Y,phi,phigx,phigy,phigxy,j,i,X0(1),X0(2),2); 
                               H(hx,hy,X,Y,phi,phigx,phigy,phigxy,j,i,X0(1),X0(2),3)];
                
            phigxn1(j,i) = phigradn1(1);
            phigyn1(j,i) = phigradn1(2);
                        
        end
    end
    
    % Update Level-Set & Gradient!
    
    phi = phin1;
    phigx = phigxn1;
    phigy = phigyn1;
    
    % Produce redfined surface for tracking contour
    if (plotstyle == 1 || plotstyle == 3)

        for j = 1:yr+1
            for i = 1:xr+1
       
                phifine(j,i) = H(hx,hy,X,Y,phi,phigx,phigy,phigxy,(idivide(int32(j-1),int32(refine))+1), ...
                            (idivide(int32(i-1),int32(refine))+1),XR(j,i),YR(j,i),1);
        
            end
        end
        
    end
    
    if (plotstyle == 5)
        for j = 1:yr+1
            for i = 1:xr+1
       
                phigxi(j,i) = H(hx,hy,X,Y,phi,phigx,phigy,phigxy,(idivide(int32(j-1),int32(refine))+1), ...
                            (idivide(int32(i-1),int32(refine))+1),XR(j,i),YR(j,i),2);
                phigyi(j,i) = H(hx,hy,X,Y,phi,phigx,phigy,phigxy,(idivide(int32(j-1),int32(refine))+1), ...
                            (idivide(int32(i-1),int32(refine))+1),XR(j,i),YR(j,i),3);
        
            end
        end
    end
    
    
    % Plot
    
    figure(1)
    
    if plotstyle == 1        % Refined Cubic Interpolation (Better Contour)
        contour(XR,YR,phifine,[0 0],'w','LineWidth',2)
        surface(XR,YR,phifine)
        view([-.7,-.7,1]);
    elseif plotstyle == 2    % View using only calculated points of phi (Poor contour)
        contour(X,Y,phi,[0 0],'w','LineWidth',2)
        surface(X,Y,phi)
        view([-.7,-.7,1]);
    elseif plotstyle == 3    % 2D with refined contour
        %quiver(x,y,phigx, phigy)
        %hold on
        contourf(XR,YR,phifine0,[0 0],'g','LineWidth',2)
        hold on
        contour(XR,YR,phifine,[0 0],'r','LineWidth',2)
    elseif plotstyle == 4
        quiver(x,y,phigx, phigy)
        hold on
        contour(x,y,phigxy)
        hold on
        contour(X,Y,phi,[0 0],'r','LineWidth',2)
    elseif plotstyle == 5
        contour(X,Y,phigxy)
        hold on
        contour(X,Y,phi,[0 0],'r','LineWidth',2)
        hold on
        quiver(XR,YR,phigxi, phigyi)
        hold on
        quiver(x,y,phigx, phigy, 'k')
    else                % 2D with poor contour
        quiver(x,y,phigx, phigy)
        hold on
        contour(X,Y,phi,[0 0],'r','LineWidth',2)
    end
    
    hold on
    axis equal
    axis([0 1 0 1])
    pic = figure(1);
    aviobjy = addframe(aviobjy,pic);
    clf(1)
    
    time
    
end

close all
aviobjy = close(aviobjy);

end











function result = H(hx,hy,X,Y,phi,phigx,phigy,phigxy,j,i,x,y,flavor)

% Input flavor: 1 = Level-set, 2 = x-derv, 3 = y-derv
% Input indices of lower-left node (of the cell) & desired coordinate

% Ensure that the bottom-left node is used in cell interpolation.
if (x < X(j,i) && y < Y(j,i))
    
    i = i - 1;
    j = j - 1;
    
elseif (x < X(j,i))
        
    i = i - 1;
    
elseif (y < Y(j,i))
    
    j = j - 1;
    
end

% Prevent interpolation using top and right boundary nodes.
if (j == size(X,1))
    
    j = j - 1;
    
end

if (i == size(X,2))
    
    i = i - 1;
    
end
    
% Interpolating variables & functions
chi = (x - X(j,i))./hx;
eta = (y - Y(j,i))./hy;

f = @(x) 1 - 3.*x.^2 + 2.*x.^3;
g = @(x) x.*(1 - x).^2;

fprime = @(x) 6.*(x.^2 - x);
gprime = @(x) 1 - 4.*x + 3.*x.^2;

% Level-Set Interpolation
if (flavor == 1)

    result = phi(j,i).*f(chi).*f(eta) + phi(j+1,i).*f(chi).*f(1-eta) ...
           + phi(j,i+1).*f(1-chi).*f(eta) + phi(j+1,i+1).*f(1-chi).*f(1-eta) ...
      + hy.*(phigy(j,i).*f(chi).*g(eta) - phigy(j+1,i).*f(chi).*g(1-eta) ...
           + phigy(j,i+1).*f(1-chi).*g(eta) - phigy(j+1,i+1).*f(1-chi).*g(1-eta)) ...
      + hx.*(phigx(j,i).*g(chi).*f(eta) + phigx(j+1,i).*g(chi).*f(1-eta) ...
           - phigx(j,i+1).*g(1-chi).*f(eta) - phigx(j+1,i+1).*g(1-chi).*f(1-eta)) ...
  + hx.*hy.*(phigxy(j,i).*g(chi).*g(eta) - phigxy(j+1,i).*g(chi).*g(1-eta) ...
           - phigxy(j,i+1).*g(1-chi).*g(eta) + phigxy(j+1,i+1).*g(1-chi).*g(1-eta));

        
         
% Level-Set Gradient Interpolation (x-direction)         
elseif (flavor == 2)
   
    result = 1./hx.*(phi(j,i).*fprime(chi).*f(eta) + phi(j+1,i).*fprime(chi).*f(1-eta) ...
                   - phi(j,i+1).*fprime(1-chi).*f(eta) - phi(j+1,i+1).*fprime(1-chi).*f(1-eta)) ...
          + hy./hx.*(phigy(j,i).*fprime(chi).*g(eta) - phigy(j+1,i).*fprime(chi).*g(1-eta) ...
                   - phigy(j,i+1).*fprime(1-chi).*g(eta) + phigy(j+1,i+1).*fprime(1-chi).*g(1-eta)) ...
                  + (phigx(j,i).*gprime(chi).*f(eta) + phigx(j+1,i).*gprime(chi).*f(1-eta) ...
                   + phigx(j,i+1).*gprime(1-chi).*f(eta) + phigx(j+1,i+1).*gprime(1-chi).*f(1-eta)) ...
              + hy.*(phigxy(j,i).*gprime(chi).*g(eta) - phigxy(j+1,i).*gprime(chi).*g(1-eta) ...
                   + phigxy(j,i+1).*gprime(1-chi).*g(eta) - phigxy(j+1,i+1).*gprime(1-chi).*g(1-eta));

                                 
% Level-set Gradient Interpolation (y-direction)
elseif (flavor == 3)
    
    result = 1./hy.*(phi(j,i).*f(chi).*fprime(eta) - phi(j+1,i).*f(chi).*fprime(1-eta) ...
                   + phi(j,i+1).*f(1-chi).*fprime(eta) - phi(j+1,i+1).*f(1-chi).*fprime(1-eta)) ...
                  + (phigy(j,i).*f(chi).*gprime(eta) + phigy(j+1,i).*f(chi).*gprime(1-eta) ...
                   + phigy(j,i+1).*f(1-chi).*gprime(eta) + phigy(j+1,i+1).*f(1-chi).*gprime(1-eta)) ...
          + hx./hy.*(phigx(j,i).*g(chi).*fprime(eta) - phigx(j+1,i).*g(chi).*fprime(1-eta) ...
                   - phigx(j,i+1).*g(1-chi).*fprime(eta) + phigx(j+1,i+1).*g(1-chi).*fprime(1-eta)) ...
              + hx.*(phigxy(j,i).*g(chi).*gprime(eta) + phigxy(j+1,i).*g(chi).*gprime(1-eta) ...
                   - phigxy(j,i+1).*g(1-chi).*gprime(eta) - phigxy(j+1,i+1).*g(1-chi).*gprime(1-eta));
                 
end

end













function result = Hi(hx,hy,X,Y,x,y,time,flavor)

% Velocity (Known ONLY at gridpoints)
T = 8;
ufun = @(x,y,time) - 2.*cos(pi.*time./T).*(sin(pi.*x).^2).*sin(pi.*y).*cos(pi.*y);
vfun = @(x,y,time)   2.*cos(pi.*time./T).*(sin(pi.*y).^2).*sin(pi.*x).*cos(pi.*x);

% Find nearest left-bottom node
xi = round((x - mod(x,hx))/(hx)) + 1;
yi = round((y - mod(y,hy))/(hy)) + 1;
% Watch out for edges
if (xi > 1/hx)% round(1/hx + 1))
    xi = xi - 1;
end
if (yi > 1/hy)% round(1/hy + 1))
    yi = yi - 1;
end


if (flavor == 1 || flavor == 4)
        
    l1 = 1/hx*[  (x - X(yi,xi+1))  ...
                -(x - X(yi,xi))];
        
    l2 = 1/hy*[  (y - Y(yi+1,xi))  ...
                -(y - Y(yi,xi))];
                        
elseif (flavor == 2 || flavor == 5)
    
    l1 = 1/hx*[  1  ...
                -1];
        
    l2 = 1/hy*[  (y - Y(yi+1,xi))  ...
                -(y - Y(yi,xi))];
        
elseif (flavor == 3 || flavor == 6)
    
    l1 = 1/hx*[  (x - X(yi,xi+1))  ...
                -(x - X(yi,xi))];
        
    l2 = 1/hy*[  1  ...
                -1];
                        
end
    
    
if (flavor < 4) % x-velocity and derivatives
    
    result = ufun(X(yi,xi),     Y(yi,xi),     time)   *l1(1)   *l2(1)    ...
           + ufun(X(yi,xi+1),   Y(yi,xi+1),   time)   *l1(2)   *l2(1)    ...
           + ufun(X(yi+1,xi),   Y(yi+1,xi),   time)   *l1(1)   *l2(2)    ...
           + ufun(X(yi+1,xi+1), Y(yi+1,xi+1), time)   *l1(2)   *l2(2);
                         
       
       
else            % y-velocity and derivatives
            
    result = vfun(X(yi,xi),     Y(yi,xi),     time)   *l1(1)   *l2(1)    ...
           + vfun(X(yi,xi+1),   Y(yi,xi+1),   time)   *l1(2)   *l2(1)    ...
           + vfun(X(yi+1,xi),   Y(yi+1,xi),   time)   *l1(1)   *l2(2)    ...
           + vfun(X(yi+1,xi+1), Y(yi+1,xi+1), time)   *l1(2)   *l2(2);
    
end

end
















function result = Hi4(hx,hy,X,Y,x,y,time,flavor)

% Velocity (Known ONLY at gridpoints)
T = 8;
ufun = @(x,y,time) - 2.*cos(pi.*time./T).*(sin(pi.*x).^2).*sin(pi.*y).*cos(pi.*y);
vfun = @(x,y,time)   2.*cos(pi.*time./T).*(sin(pi.*y).^2).*sin(pi.*x).*cos(pi.*x);

% Find nearest left-bottom node
xi = round((x - mod(x,hx))/(hx)) + 1;
yi = round((y - mod(y,hy))/(hy)) + 1;
% Find nearest interpolating left-bottom node (using 4x4 stencil)
if (xi >= 1/hx)% round(1/hx + 1))
    xi = xi - 2;
end
if (yi >= 1/hy)% round(1/hy + 1))
    yi = yi - 2;
end
xi = round(xi - mod(xi-1,3));
yi = round(yi - mod(yi-1,3));


if (flavor == 1 || flavor == 4)
        
    l1 = 1/(6*hx^3)*[  -(x - X(yi,xi+1))*(x - X(yi,xi+2))  *(x - X(yi,xi+3))    ...
                      3*(x - X(yi,xi))  *(x - X(yi,xi+2))  *(x - X(yi,xi+3))    ...
                     -3*(x - X(yi,xi))  *(x - X(yi,xi+1))  *(x - X(yi,xi+3))    ...
                        (x - X(yi,xi))  *(x - X(yi,xi+1))  *(x - X(yi,xi+2))];
        
    l2 = 1/(6*hy^3)*[  -(y - Y(yi+1,xi))*(y - Y(yi+2,xi))  *(y - Y(yi+3,xi))    ...
                      3*(y - Y(yi,xi))  *(y - Y(yi+2,xi))  *(y - Y(yi+3,xi))    ...
                     -3*(y - Y(yi,xi))  *(y - Y(yi+1,xi))  *(y - Y(yi+3,xi))    ...
                        (y - Y(yi,xi))  *(y - Y(yi+1,xi))  *(y - Y(yi+2,xi))];
                        
elseif (flavor == 2 || flavor == 5)
    
    l1 = 1/(6*hx^3)*[  -((x - X(yi,xi+1))*(x - X(yi,xi+2)) + (x - X(yi,xi+2))*(x - X(yi,xi+3)) + (x - X(yi,xi+1))*(x - X(yi,xi+3)))    ...
                      3*((x - X(yi,xi))*(x - X(yi,xi+2))   + (x - X(yi,xi+2))*(x - X(yi,xi+3)) + (x - X(yi,xi))*(x - X(yi,xi+3)))    ...
                     -3*((x - X(yi,xi))*(x - X(yi,xi+1))   + (x - X(yi,xi+1))*(x - X(yi,xi+3)) + (x - X(yi,xi))*(x - X(yi,xi+3)))    ...
                        ((x - X(yi,xi))*(x - X(yi,xi+1))   + (x - X(yi,xi+1))*(x - X(yi,xi+2)) + (x - X(yi,xi))*(x - X(yi,xi+2)))];
        
    l2 = 1/(6*hy^3)*[  -(y - Y(yi+1,xi))*(y - Y(yi+2,xi))  *(y - Y(yi+3,xi))    ...
                      3*(y - Y(yi,xi))  *(y - Y(yi+2,xi))  *(y - Y(yi+3,xi))    ...
                     -3*(y - Y(yi,xi))  *(y - Y(yi+1,xi))  *(y - Y(yi+3,xi))    ...
                        (y - Y(yi,xi))  *(y - Y(yi+1,xi))  *(y - Y(yi+2,xi))];
        
elseif (flavor == 3 || flavor == 6)
    
    l1 = 1/(6*hx^3)*[  -(x - X(yi,xi+1))*(x - X(yi,xi+2))  *(x - X(yi,xi+3))    ...
                      3*(x - X(yi,xi))  *(x - X(yi,xi+2))  *(x - X(yi,xi+3))    ...
                     -3*(x - X(yi,xi))  *(x - X(yi,xi+1))  *(x - X(yi,xi+3))    ...
                        (x - X(yi,xi))  *(x - X(yi,xi+1))  *(x - X(yi,xi+2))];

    l2 = 1/(6*hy^3)*[  -((y - Y(yi+1,xi))*(y - Y(yi+2,xi)) + (y - Y(yi+2,xi))*(y - Y(yi+3,xi)) + (y - Y(yi+1,xi))*(y - Y(yi+3,xi)))    ...
                      3*((y - Y(yi,xi))  *(y - Y(yi+2,xi)) + (y - Y(yi+2,xi))*(y - Y(yi+3,xi)) + (y - Y(yi,xi))*(y - Y(yi+3,xi)))    ...
                     -3*((y - Y(yi,xi))  *(y - Y(yi+1,xi)) + (y - Y(yi+1,xi))*(y - Y(yi+3,xi)) + (y - Y(yi,xi))*(y - Y(yi+3,xi)))    ...
                        ((y - Y(yi,xi))  *(y - Y(yi+1,xi)) + (y - Y(yi+1,xi))*(y - Y(yi+2,xi)) + (y - Y(yi,xi))*(y - Y(yi+2,xi)))];
                        
end
    
    
if (flavor < 4) % x-velocity and derivatives
    
    result = ufun(X(yi,xi),    Y(yi,xi),    time)   *l1(1)    *l2(1)    ...
           + ufun(X(yi,xi+1),  Y(yi,xi+1),  time)   *l1(2)    *l2(1)    ...
           + ufun(X(yi,xi+2),  Y(yi,xi+2),  time)   *l1(3)    *l2(1)    ...
           + ufun(X(yi,xi+3),  Y(yi,xi+3),  time)   *l1(4)    *l2(1)    ...
           + ufun(X(yi+1,xi),  Y(yi+1,xi),  time)   *l1(1)    *l2(2)    ...
           + ufun(X(yi+1,xi+1),Y(yi+1,xi+1),time)   *l1(2)    *l2(2)    ...
           + ufun(X(yi+1,xi+2),Y(yi+1,xi+2),time)   *l1(3)    *l2(2)    ...
           + ufun(X(yi+1,xi+3),Y(yi+1,xi+3),time)   *l1(4)    *l2(2)    ...
           + ufun(X(yi+2,xi),  Y(yi+2,xi),  time)   *l1(1)    *l2(3)    ...
           + ufun(X(yi+2,xi+1),Y(yi+2,xi+1),time)   *l1(2)    *l2(3)    ...
           + ufun(X(yi+2,xi+2),Y(yi+2,xi+2),time)   *l1(3)    *l2(3)    ...
           + ufun(X(yi+2,xi+3),Y(yi+2,xi+3),time)   *l1(4)    *l2(3)    ...
           + ufun(X(yi+3,xi),  Y(yi+3,xi),  time)   *l1(1)    *l2(4)    ...
           + ufun(X(yi+3,xi+1),Y(yi+3,xi+1),time)   *l1(2)    *l2(4)    ...
           + ufun(X(yi+3,xi+2),Y(yi+3,xi+2),time)   *l1(3)    *l2(4)    ...
           + ufun(X(yi+3,xi+3),Y(yi+3,xi+3),time)   *l1(4)    *l2(4);
                         
       
       
else            % y-velocity and derivatives
            
    result = vfun(X(yi,xi),    Y(yi,xi),    time)   *l1(1)    *l2(1)    ...
           + vfun(X(yi,xi+1),  Y(yi,xi+1),  time)   *l1(2)    *l2(1)    ...
           + vfun(X(yi,xi+2),  Y(yi,xi+2),  time)   *l1(3)    *l2(1)    ...
           + vfun(X(yi,xi+3),  Y(yi,xi+3),  time)   *l1(4)    *l2(1)    ...
           + vfun(X(yi+1,xi),  Y(yi+1,xi),  time)   *l1(1)    *l2(2)    ...
           + vfun(X(yi+1,xi+1),Y(yi+1,xi+1),time)   *l1(2)    *l2(2)    ...
           + vfun(X(yi+1,xi+2),Y(yi+1,xi+2),time)   *l1(3)    *l2(2)    ...
           + vfun(X(yi+1,xi+3),Y(yi+1,xi+3),time)   *l1(4)    *l2(2)    ...
           + vfun(X(yi+2,xi),  Y(yi+2,xi),  time)   *l1(1)    *l2(3)    ...
           + vfun(X(yi+2,xi+1),Y(yi+2,xi+1),time)   *l1(2)    *l2(3)    ...
           + vfun(X(yi+2,xi+2),Y(yi+2,xi+2),time)   *l1(3)    *l2(3)    ...
           + vfun(X(yi+2,xi+3),Y(yi+2,xi+3),time)   *l1(4)    *l2(3)    ...
           + vfun(X(yi+3,xi),  Y(yi+3,xi),  time)   *l1(1)    *l2(4)    ...
           + vfun(X(yi+3,xi+1),Y(yi+3,xi+1),time)   *l1(2)    *l2(4)    ...
           + vfun(X(yi+3,xi+2),Y(yi+3,xi+2),time)   *l1(3)    *l2(4)    ...
           + vfun(X(yi+3,xi+3),Y(yi+3,xi+3),time)   *l1(4)    *l2(4);
    
end

end








