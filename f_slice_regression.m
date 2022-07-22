function [x, y, r, theta, out_counter, yy]=f_slice_regression(x0, y0, r0, theta0, r_rocket, h, rr0, fps, p, A_t, A0, pc0, rho_b, c_star, a, n)

delta_theta=pi/p; %[rad] angle occupied by a star point "slice" (slice=half star point)
R=[0 1; -1 0]; %90 deg clockwise rotation matrix
x(1,:)=x0; %x coordinates vector initialization
y(1,:)=y0; %y coordinates vector initialization
r(1,:)=r0; %radial coordinates vector initialization
theta(1,:)=theta0; %angular coordinates vector initialization

yy(1,:)=0; %yy coordinates (=radial combustion coordinate) vector initialization (for analytical method)

A_b(1,:)=A0; %combustion surfaces vector initialization
pc(1,:)=pc0; %chamber pressures vector initialization
rr=rr0; %regression ratio initialization

stop_counter=zeros(1,size(r,2)); %control vector for nodes outside combustion chamber
out_counter=zeros(1,size(r,2)); %control vector for nodes outside star point slice perimeter

i=2; % IMPORTANT: "i" is the temporal iteration index
while ~isequal(stop_counter,ones(1,size(r,2)))
   
   %iterative cycle for nodes regression
   for j=1:1:size(r,2)
       
       %IMPORTANT: "j" is the index indicating which node we are delaing with
       
       if j==1 %special case: node 1 
           r(i,j)=r(i-1,j)+(rr/fps); %increasing radial coordinate by rr/fps
           theta(i,j)=theta(i-1,j); %theta is constant for first node
           x(i,j)=r(i,j)*cos(theta(i-1,j)); %calculating x coordinate
           y(i,j)=r(i,j)*sin(theta(i-1,j)); %calculating y coordinate
           if stop_counter(j)==1 %checking that the node wasn't already outside comb chamb
               %if it was out already, then its coordinates must not change
               r(i,j)=r(i-1,j);
               theta(i,j)=theta(i-1,j);
               x(i,j)=x(i-1,j);
               y(i,j)=y(i-1,j);
           elseif r(i,j)>=r_rocket %checking that the node did't exit comb chamb at current iteration
               %if it did, then I move it on the comb chamb perimeter
               x(i,j)=r_rocket*cos(theta(i-1,j));
               y(i,j)=r_rocket*sin(theta(i-1,j));
               stop_counter(j)=1;
           end
           
       elseif j==size(r,2) %special case: last node   
           vec=[x(i-1,j)-x(i-1,j-1) ; y(i-1,j)-y(i-1,j-1)]; %calculating tangent vector to node
           vec_rot=R*vec; %rotating tangent vec to obtain normal vec
           multiplier=rr/(fps*norm(vec_rot));
           %increasing x and y coordinates
           x(i,j)=x(i-1,j)+(vec_rot(1,1)*multiplier); 
           y(i,j)=y(i-1,j)+(vec_rot(2,1)*multiplier); 
           coordinate=[x(i,j) y(i,j)];
           r(i,j)=norm(coordinate); %calculating new radial coordinate
           %calculating new angular coordinate
           if x(i,j)<0
                theta(i,j)=atan(y(i,j)/x(i,j))+pi;
           else
                theta(i,j)=atan(y(i,j)/x(i,j));
           end
           if stop_counter(j)==1 %checking that the node wasn't already outside comb chamb
               %if it was out already, then its coordinates must not change
               r(i,j)=r(i-1,j);
               theta(i,j)=theta(i-1,j);
               x(i,j)=x(i-1,j);
               y(i,j)=y(i-1,j);
           elseif r(i,j)>=r_rocket %checking that the node did't exit comb chamb at current iteration
               %if it did, then I move it on the comb chamb perimeter
               x(i,j)=r_rocket*cos(theta(i-1,j));
               y(i,j)=r_rocket*sin(theta(i-1,j));
               stop_counter(j)=1;
           end
           if theta(i,j)>delta_theta %checking that the node didn't exit star point slice perimeter at current iteration
               out_counter(i,j)=1;
               stop_counter(j)=1;
           end
                  
       else %for all other nodes (sames as previous elseif block)
           vec=[x(i-1,j+1)-x(i-1,j-1) ; y(i-1,j+1)-y(i-1,j-1)];
           vec_rot=R*vec;
           multiplier=rr/(fps*norm(vec_rot));
           x(i,j)=x(i-1,j)+(vec_rot(1,1)*multiplier);
           y(i,j)=y(i-1,j)+(vec_rot(2,1)*multiplier);
           coordinate=[x(i,j) y(i,j)];
           r(i,j)=norm(coordinate);
           if x(i,j)<0
                theta(i,j)=atan(y(i,j)/x(i,j))+pi;
           else
                theta(i,j)=atan(y(i,j)/x(i,j));
           end
           if stop_counter(j)==1
               r(i,j)=r(i-1,j);
               theta(i,j)=theta(i-1,j);
               x(i,j)=x(i-1,j);
               y(i,j)=y(i-1,j);
           elseif r(i,j)>= r_rocket
               x(i,j)=r_rocket*cos(theta(i-1,j));
               y(i,j)=r_rocket*sin(theta(i-1,j));
               stop_counter(j)=1;
           end
           if theta(i,j)>delta_theta
               out_counter(i,j)=1;
               stop_counter(j)=1;
           end        
       end
   end
   
   %calculating radial combustion coordinate at current i-iteration
   yy(1,i)=yy(1,i-1)+(rr/fps);
   
   %calculating current A_b and regression ratio for next i-iteration
   S_b=0; %perimeter initialization
   in_vec=find(out_counter(i,:)==0); %checking which nodes to consider
   x_count=x(i,in_vec); %extracting x coordinates
   y_count=y(i,in_vec); %extracting y coordinates
   for j=2:size(x_count,2) %calculating combustion perimeter at current i-iteration
       if r(i,j)< r_rocket
           S_b=S_b+sqrt(((x_count(j)-x_count(j-1))^2)+((y_count(j)-y_count(j-1))^2));
       end
   end
   A_b(1,i)=S_b*2*p*h; %calculating combustion surface from perimeter
   K=A_b(1,i)/A_t; %new Klemmung factor
   pc(1,i)=(K*rho_b*c_star*a)^(1/(1-n)); % [Pa] new chamber pressure   
   rr=a*(pc(1,i)^n); % [m/s] regression ratio for next iteration

   
i=i+1;
end

end