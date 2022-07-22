function [x0, y0, r0, theta0, A0, pc0, rr0]=f_initial_slice(n_tip, n_fillet, n_side, p, eps, r_tip, r_fillet, z, h, A_t, rho_b, c_star, a, n) 

%NOTE: for details on geometry terminology see attached scheme

%TIP NODES
delta_theta=pi/p; %[rad] angle occupied by a star point "slice" (slice=half star point)
delta_theta_side=delta_theta*eps; %[rad] angle occupied by star point slice side
delta_theta_tip=delta_theta-delta_theta_side; %[rad] angle occupied by star point slice tip

r0_tip=r_tip*ones(1, n_tip); %radial coordinates of star point slice tip nodes
theta0_tip=linspace(0,delta_theta_tip, n_tip); %angular coordinates of star point slice tip nodes
 
%FILLET NODES
th=z*delta_theta; %star point throat half-opening (th="theta half")
xC_fillet=(r0_tip(end)-r_fillet)*cos(theta0_tip(end)); %fillet centre x coordinate
yC_fillet=(r0_tip(end)-r_fillet)*sin(theta0_tip(end)); %fillet centre y coordinate
thetaC_fillet_end=(pi/2)+delta_theta_side-th; %fillet angle with respect to its centre
delta_thetaC_fillet=thetaC_fillet_end/n_fillet; %angle between fillet nodes with respect to its centre
thetaC_fillet=linspace(theta0_tip(end)+delta_thetaC_fillet,theta0_tip(end)+thetaC_fillet_end, n_fillet); %angular coordinates of fillet nodes with respect to its centre

for i=1:n_fillet %fillet nodes coordinates with respect to grain centre
   x0_fillet(i)=r_fillet*cos(thetaC_fillet(i))+xC_fillet;
   y0_fillet(i)=r_fillet*sin(thetaC_fillet(i))+yC_fillet;
   r0_fillet(i)=sqrt((x0_fillet(i)^2)+(y0_fillet(i)^2));
   theta0_fillet(i)=atan(y0_fillet(i)/x0_fillet(i));
end

%STAR POINT SIDE NODES
gamma=th-(delta_theta-theta0_fillet(end));
r_base=r0_fillet(end)*(sin(gamma)/sin(pi-th)); %throat radius calculated with sine-law
x_side_start=r0_fillet(end)*cos(theta0_fillet(end)); %x coordinate of star point side beginning (from fillet)
y_side_start=r0_fillet(end)*sin(theta0_fillet(end)); %y coordinate of star point side beginning (from fillet)
x_side_end=r_base*cos(delta_theta); %x coordinate of star point side end (at throat point)
y_side_end=r_base*sin(delta_theta); %y coordinate of star point side end (at throat point)

dx=(x_side_end-x_side_start)/n_side; %x-wise distance between nodes on star point side
dy=(y_side_end-y_side_start)/n_side; %y-wise distance between nodes on star point side

for i=1:n_side %coordinates for all nodes on star point side
    x_side(i)=x_side_start+(i*dx);
    y_side(i)=y_side_start+(i*dy);
    if x_side(i)<0
        theta0_side(i)=atan(y_side(i)/x_side(i))+pi;
    else
        theta0_side(i)=atan(y_side(i)/x_side(i));
    end
    r0_side(i)=sqrt((x_side(i)^2)+(y_side(i)^2));
end

%creating single vectors with all angles and all radii
r0=[r0_tip r0_fillet r0_side];
theta0=[theta0_tip theta0_fillet theta0_side];
x0=r0.*cos(theta0);
y0=r0.*sin(theta0);

%initial combustion surface A0
S0=0; %perimeter initialization
for j=2:size(x0,2) %calculating perimeter
    S0=S0+sqrt(((x0(j)-x0(j-1))^2)+((y0(j)-y0(j-1))^2));
end
A0=S0*2*p*h; %initial combustion surface

%initial combustion parameters
K0=A0/A_t;                          % initial Klemmung factor
pc0=(K0*rho_b*c_star*a)^(1/(1-n));  % [Pa] initial chamber pressure
rr0=a*(pc0^n);                      % [m/s] initial regression ratio

end