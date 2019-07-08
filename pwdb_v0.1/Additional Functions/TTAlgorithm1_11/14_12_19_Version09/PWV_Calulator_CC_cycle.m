function t_int = PWV_Calulator_CC_cycle(t,signal,indmaxs,indmins,gradins)

dt = t(2)-t(1);

% Determine most approximate length of cycle, most consistent kernel, and
% the approximate transit time
dCycle1 = abs(diff(gradins(1,:,1)));
dCycle2 = abs(diff(gradins(2,:,1)));

diff1 = std(dCycle1);
diff2 = std(dCycle2);
if diff1 <= diff2
    CApprox = round( median(dCycle1) );
    kernel = gradins(1,:,1);
else
    CApprox = round( median(dCycle2) );
    kernel = gradins(2,:,1);
end

kernelExtend = round(0.2*CApprox);
trash = find(kernel > kernelExtend+1);
kernel = kernel(trash);

TTApprox = median(diff(gradins(:,:,1)));

if size(gradins,2) > 1
    
    t_int = [];
    
    for i=2:length(kernel)
        
        kernel1 = [kernel(i) , kernel(i-1)];
        
        M1 = min(signal(1,kernel1(1) : kernel1(2)));
        M2 = min(signal(2,kernel1(1) : kernel1(2)));
        
        sig1 = [M1*ones(1,2*TTApprox), signal(1,kernel1(1)-kernelExtend : kernel1(2)), M1*ones(1,2*TTApprox)];
        sig2 = [M2*ones(1,2*TTApprox), signal(2,kernel1(1)-kernelExtend : kernel1(2)), M2*ones(1,2*TTApprox)];
        
        % Normalise Carotid Wave Data
        trash = sig1 - min(sig1);
        sig1 = trash/max(trash);
        
        % Normalise Femoral Wave Data
        trash = sig2 - min(sig2);
        sig2 = trash/max(trash);
        
        y = [];
        x = [];
        
        for j= -round(TTApprox/4): 1 : round(1.5*TTApprox)
            
            sig1_new = sig1( (2*TTApprox+1)   : (2*TTApprox + kernelExtend + round(0.8*CApprox))   );
            sig2_new = sig2( (2*TTApprox+1)+j : (2*TTApprox + kernelExtend + round(0.8*CApprox))+j );
           
            trash = find(sig2_new);
            
            xx = sum(sig1_new.*sig2_new);
            y = [y,xx];
            x = [x,(j*dt)];
            
            if xx > max(y(1:end-1))
                max_ind = j;
                max_CC = xx;
            else
            end
        end
        
        loop=1;
        k=1;
        while loop==1
            k=k+1;
            poly = polyfit(x(k-1:k+1),y(k-1:k+1),1);
            if poly(1) > 10
                ind1 = k;
                loop = 0;
            end
            if (k+1) == length(y)
                ind1 = 1;
                loop = 0;
            end
        end
        
        [~,I] = max(y(ind1:end));
        I = I + (ind1-1);
        interp_span = 5;
        
        if length(y) < (I+interp_span)
            interp_span = length(y) - I;
        elseif (I-interp_span) < 1
            interp_span = I-1;
        end
        
        interp_int = [I-interp_span:I+interp_span];
        int = find(interp_int);
        interp_int = interp_int(int);
        
        poly = polyfit(x(interp_int),y(interp_int),2);
        x1 = [x(1):0.001:x(end)];
        y1 = poly(1).*x1.^2 + poly(2).*x1 + poly(3);
        
        m_TT = -poly(2)/(2*poly(1));
        
        t_int = [t_int, m_TT];
    end
    
else
    ind = indmins(:,:,1);
    kernel1 = [indmins(1,:,1) , indmaxs(1)];
    kernel2 = [indmins(2,:,1) , indmaxs(2)];
    
    shift = max( round( abs(kernel2(2) - kernel1(2))*4 ) , 25 );
    
    trash = round(size(signal,2)/4);
    [~,I] = min( signal(1, (kernel1(1)+trash) :end));
    max_length = I  + kernel1(1)+trash - 1;
    
    sig1_new = signal(1,:);
    
    sig2 = zeros(1,size(signal,2));
    sig2(kernel2(1):max_length) = signal(2,kernel2(1):max_length);
    
    y = [];
    x = [];
    
    for j= -round(shift):1:round(shift)
        if j<0
            sig2_new = [zeros(1,-j),sig2(1 : end-(-j)) ];
        else
            sig2_new = [sig2(j+1:end),zeros(1,j)];
        end
        
        trash = find(sig2_new);
        if numel(trash)<1
            save('error_data','')
        elseif (trash(end)-trash(1)) >= (length( (kernel2(1)+1) :max_length)-1)
            xx = sum(sig1_new.*sig2_new);
            y = [y,xx];
            x = [x,(j*dt)];
            
            if xx > max(y(1:end-1))
                max_ind = j;
                max_CC = xx;
            else
            end
        else
        end
    end
    
    % Recreate max correlation to find the Correlation coefficient
    j = max_ind;
    if j<0
        sig2_new = [zeros(1,-j),sig2(1 : end-(-j)) ];
    else
        sig2_new = [sig2(j+1:end),zeros(1,j)];
    end
    
    ind_x = find(sig2_new);
    SumProx1 = sum(sig1_new(ind_x).^2);
    SumProx2 = sum(sig2_new(ind_x).^2);
    CC_ref = max(SumProx1,SumProx2);
    
    CrossCoeff = max_CC/CC_ref;
    
    loop=1;
    i=1;
    while loop==1
        i=i+1;
        poly = polyfit(x(i-1:i+1),y(i-1:i+1),1);
        if poly(1) > 10
            ind = i;
            loop = 0;
        end
        if (i+1) == length(y)
            ind = 1;
            loop = 0;
        end
    end
    
    [~,I] = max(y(ind:end));
    I = I + (ind-1);
    interp_span = 5;
    
    if length(y) < (I+interp_span)
        interp_span = length(y) - I;
    elseif (I-interp_span) < 1
        interp_span = I-1;
    end
    
    interp_int = [I-interp_span:I+interp_span];
    int = find(interp_int);
    interp_int = interp_int(int);
    
    poly = polyfit(x(interp_int),y(interp_int),2);
    x1 = [x(1):0.001:x(end)];
    y1 = poly(1).*x1.^2 + poly(2).*x1 + poly(3);
    
    t_int = -poly(2)/(2*poly(1));
    
end