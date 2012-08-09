% Script to collect & organize specific data

start = 1;
filenum = start;
numframes = 850;
result = zeros(numframes,2);

for j = start:numframes
        
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

    % What are we looking for?  Goes here:
    
    xl = 240;
    k = 1;
    
    PHI_for_plot;
    
    while PHI_for_plot(k,xl) > 0
        k = k+1;
    end
    k = k-1;
    
    result(j,1) = j*.02;
    result(j,2) = Y_for_plot(k,xl) - PHI_for_plot(k,xl)/PHIY_for_plot(k,xl);
    
    fclose(fid);
    clf
    
    filenum = filenum +1  

end

plot(result(:,2))

