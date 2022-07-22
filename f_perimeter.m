function [x_slice, y_slice]=f_perimeter(p, r_rocket) 

delta_theta=pi/p; %[rad] angle occupied by a star point "slice" (slice=half star point)
theta_slice=linspace(0,delta_theta,100); %angular coordinates for combustion chamber perimeter nodes
x_slice=r_rocket.*cos(theta_slice); %x coordinates
y_slice=r_rocket.*sin(theta_slice); %y coordinates

end