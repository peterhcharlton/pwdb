function t_int= PWV_Calulator_Least_Squares(t,signal,indmaxs,indmins)

%%%%%%%%%%%%% Least Squares - Start %%%%%%%%%%%%%

dt = t(2)-t(1);

trash = diff(indmaxs);
trash = median(trash);

ind = find(indmins(1,:,1));
kernel1_tot = [indmins(1, : ,1);indmaxs(1, : )];
ind = find(indmins(2,:,1));
kernel2_tot = [indmins(2, : ,1);indmaxs(2, : )];
Prox = 4;

pixPCycle = round(abs(mean(diff(kernel1_tot(1,:)))));

padding1 = round(pixPCycle/3);
padding2 = round(2*pixPCycle/3);

count = 1;

if length(ind) > 1 % if there are more than one cycles
    
    for k = 2:length(ind)
        
        kernel1 = kernel1_tot(:,k);
        kernel2 = kernel2_tot(:,k);
        
        if kernel1(1)<kernel1(2) && kernel2(1)<kernel2(2)
            if kernel1(1) > kernel2(1)
                kernel1(1) = kernel2(1);
            end
            shift = round( abs(kernel1(2) - kernel1(1)));
            
            max_lengthBack = max( [(kernel1(2) - kernel1(1)) , (kernel2(2) - kernel2(1)) , abs(kernel2(2) - kernel1(1))] );
            max_lengthForward = round(max_lengthBack/2);
            
            if size(signal,2) < kernel2(2)+max_lengthForward
                max_lengthForward = round(size(signal,2) - kernel2(2));
            end
            
            % Normalise the signals and lift off the x axis
            min1 = min(signal(1,kernel1(1):kernel1(2)));
            trash1 = signal(1,:) - min1;
            sig1_norm = trash1/(trash1(kernel1(2))) + 0.05;
            
            min2 = min(signal(2,kernel2(1):kernel2(2)));
            trash1 = signal(2,:) - min2;
            sig2_norm = trash1/(trash1(kernel2(2)));
            
            sig2_base = zeros(1,length(t));
            lim1 = kernel2(2)-max_lengthBack;
            if lim1 < 1
                lim1 = 1;
            end
            lim2 = kernel2(2)+max_lengthForward;
            if lim2 > length(sig2_base)
                lim2 = length(sig2_base);
            end
            sig2_base(lim1:lim2) = sig2_norm(lim1:lim2) + 0.05;
            
            if (kernel1(1) - padding1) < 1
                padding1 = kernel1(1) - 1;
            elseif (kernel1(2) + padding2) >= length(t)
                padding2 = length(t) - kernel1(2);
            end
            
            % Signals for least squares correlation
            sig1_new = sig1_norm(kernel1(1)-padding1 : kernel1(2)+padding2);
            sig2     = sig2_base(kernel1(1)-padding1 : kernel1(2)+padding2);
            
            CCLength = length(find(sig2));
            
            y = [];
            x = [];
            
            loop = 1;
            j = -round(shift);
            
            while loop == 1
                if j<0
                    sig2_new = [zeros(1,-j),sig2(1 : end-(-j)) ];
                else
                    sig2_new = [sig2(j+1:end),zeros(1,j)];
                end
                
                trash = find(sig2_new);
                if (trash(end)-trash(1)) >= round(CCLength*0.7)
                    trash = (sig1_new.*sig2_new);
                    [C,I] = find(trash);
                    LS = 0;
                    if numel(I)>0
                        LS = sum( (sig1_new(I) -sig2_new(I)).^2 );
                    else
                        LS = 20;
                    end
                    
                    y = [y,LS];
                    x = [x,(j*dt)];
                    
                    j = j+1;
                    ind = find(sig2_new);
                    if ind < (max_lengthBack+1)
                        loop = 0;
                    end
                else
                    if isempty(y)
                        j = j+1;
                    else
                        loop = 0;
                        skip = 1;
                    end
                end
            end
            
            intZero = find(x==0);
            
            if ~isempty(intZero)
                
                trash1 = diff(y);
                MaxMinInd = find(trash1(1:end-1).*trash1(2:end)<0)+1;
                if length(MaxMinInd)>1
                    MaxMin = y(MaxMinInd);
                    startGrad = diff(MaxMin(1:2));
                    if startGrad > 0
                        trash2 = MaxMinInd(1:2:end);
                    else
                        trash2 = MaxMinInd(2:2:end);
                    end
                else
                    trash2 = MaxMinInd;
                end
                
                trash3 = abs( (intZero+Prox - trash2));
                trash4 = (trash3) .* y(trash2);
                [~,ind] = min(trash4);
                I = trash2(ind);
                
                interp_span = 1;
                if (I+interp_span) > length(y)
                    I = length(y) - interp_span;
                elseif (I-interp_span) < 1
                    I = 1 + interp_span;
                end
                interp_int = [I-interp_span:I+interp_span];
                int = find(interp_int>0);
                interp_int = interp_int(int);
                
                poly = polyfit(x(interp_int),y(interp_int),2);
                x1 = [x(1):0.001:x(end)];
                y1 = poly(1).*x1.^2 + poly(2).*x1 + poly(3);
                
                m_TT = -poly(2)/(2*poly(1));
                
                t_int(count) = m_TT;
                count = count+1;
                skip = 0;
            else
                t_int(count) = nan;
                count = count+1;
            end
        end
    end
    
else
    
    ind = indmins(:,:,1);
    kernel1 = [indmins(1,:,1) , indmaxs(1)];
    kernel2 = [indmins(2,:,1) , indmaxs(2)];
    
    shift = max( round( abs(kernel1(2) - kernel1(2))*5 ),20 );
    
    max_length = max( [(kernel1(2) - kernel1(1)) , (kernel2(2) - kernel2(1))] );
    
    trash1 = signal(1,:) - signal(1,kernel1(1));
    sig1_norm = trash1/(trash1(kernel1(2))) + 0.05;
    
    trash1 = signal(2,:) - signal(2,kernel2(1));
    sig2_norm = trash1/(trash1(kernel2(2)));
    
    sig2_base = zeros(1,length(t));
    sig2_base(kernel2(2)-max_length:kernel2(2)) = sig2_norm(kernel2(2)-max_length:kernel2(2)) + 0.05;
    
    y = [];
    x = [];
    
    loop = 1;
    j = -round(shift/4);
    while loop == 1
        if j<0
            sig2_new = [zeros(1,-j),sig2_base(1 : end-(-j)) ];
        else
            sig2_new = [sig2_base(j+1:end),zeros(1,j)];
        end
        
        trash = find(sig2_new);
        if numel(trash)<1
            loop = 0;
        elseif (trash(end)-trash(1)) >= (kernel2(2) - kernel2(1)-1)
            trash = (sig1_norm.*sig2_new);
            [C,I] = find(trash);
            LS = 0;
            if numel(I)>0
                LS = sum( (sig1_norm(I) -sig2_new(I)).^2 );
            else
                LS = 20;
            end
            
            y = [y,LS];
            x = [x,(j*dt)];
            
            j = j+1;
            ind = find(sig2_new);
            if ind < (max_length+1)
                loop = 0;
            end
        else
        end
    end
    
    [~,I] = min(y);
    interp_span = 5;
    if (I+interp_span) > length(y)
        I = length(y) - interp_span;
    elseif (I-interp_span) < 1
        I = 1 + interp_span;
    end
    interp_int = [I-interp_span:I+interp_span];
    int = find(interp_int>0);
    interp_int = interp_int(int);
    
    poly = polyfit(x(interp_int),y(interp_int),2);
    
    t_int = -poly(2)/(2*poly(1));
end