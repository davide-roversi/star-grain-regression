function [A_b_an, yy_an, phase_B, phase_slivers, iteration_B, iteration_slivers]=f_analytical_method(loops, p, eps, th, w, f, l, h, A_t, rr0, pc0, a, n, rho_b, c_star, fps)

rr=rr0;
pc=pc0;
yy_an(1,:)=0;
const=sin((eps*pi)/p); %useful variable for following operations
A_b_an(1)=2*p*l*h*((const/sin(th))+((yy_an(1,1)+f)/l)*((pi/2)+(pi/p)-th-cot(th))+((1-eps)*(pi/p)));

phase_B=0;
phase_slivers=0;

%calculating combustion surfaces with analytical formulas
for i=2:loops
    
    K=A_b_an(i-1)/A_t; % previous Klemmung factor
    pc=(K*rho_b*c_star*a)^(1/(1-n)); % [Pa] chamber pressure
    rr=a*(pc^n); % [m/s] current iteration regression ratio
    yy_an(1,i)=yy_an(1,i-1)+(rr/fps); %new radial combustion coordinate
    
    if ((yy_an(1,i)+f)/l)<=(const/cos(th))
        %phase A formula
        A_b_an(i)=2*p*l*h*((const/sin(th))+((yy_an(1,i)+f)/l)*((pi/2)+(pi/p)-th-cot(th))+((1-eps)*(pi/p)));
    else
        %phase B formula
        A_b_an(i)=2*p*l*h*((((yy_an(1,i)+f)/l)*((pi/p)+asin((l*const)/(yy_an(1,i)+f))))+((1-eps)*(pi/p)));
        if phase_B==0
           phase_B=yy_an(1,i); %saving phase B transition coordinate
           iteration_B=i;
        end
        if yy_an(1,i)>w
            if phase_slivers==0
                phase_slivers=yy_an(1,i); %saving slivers transition coordinate
                iteration_slivers=i;
            end
        end
    end
    
end

end
