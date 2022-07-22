%STAR GRAIN CONFIGURATION

clear all
clc

%-INPUT GRAIN GEOMETRY----------------------------------------------------
n_tip=50;       % number of nodes on star tip
n_fillet=200;   % number of nodes on star tip fillet
n_side=200;     % number of nodes on star side
p=8;            % number of star points
eps=0.99;       % side angle to tip angle ratio
th=37;          % [deg] half-opening of star points throat
f=0.030;        % [m] fillet radius
l=0.920;        % [m] star point length
w=0.365;        % [m] web thickness
h=6.100;        % [m] combustion chamber height
fps=0.1;        % number of surface updates per second (true animation fps)
%-------------------------------------------------------------------------

%-INPUT PROPELLANT DATA==-------------------------------------------------
R=8.3145;       % [J/mol*K] perfect gases universal constant
gamma=1.26;     % adiabatic expansion constant for exhaust gases
rho_b=1740;     % [kg/m^3] ambient temperature density of propellant
MM=0.024;       % [Kg/mol] molar mass of exhaust gases
Tc=3000;        % [K] combustion chamber temperarture
n=0.304;        % combustion coefficient (for rr in [m/s] and pc in [Pa])
a=5.170e-05;    % temperature coefficient (for rr in [m/s] and pc in [Pa])
%-------------------------------------------------------------------------

%-INPUT NOZZLE DATA-------------------------------------------------------
r_t=0.190;      % [m] nozzle throat section radius
r_e=0.860;      % [m] nozzle exit section radius
p_amb=101325;   % [Pa] ambient pressure
%-------------------------------------------------------------------------

tic

%------------------------------CALCULATIONS------------------------------%
%initializing useful variables
th=th*(pi/180);             % [rad] half-opening of star points throat
r_fillet=f;                 % [m] fillet radius
r_tip=l+f;                  % [m] star point tip radius
r_rocket=l+f+w;             % [m] combustion chamber radius
delta_theta=pi/p;           % [rad] angle occupied by a star point "slice" (slice=half star point)
z=th/delta_theta;           % star point throat half-opening to delta_theta ratio
Gam=sqrt(gamma*((2/(gamma+1))^((gamma+1)/(gamma-1))));
c_star=sqrt(Tc*R/MM)/Gam;   % [m/s] characteristic velocity
A_t=pi*(r_t)^2;             % [m^2] nozzle throat section area
A_e=pi*(r_e)^2;             % [m^2] nozzle exit section area
E=A_e/A_t;                  % area ratio

%exit Mach and pressure ratio
j=1;
fcn = @(M_e) (1./M_e)*((2/(gamma+1)).*(1+(((gamma-1)/2)*M_e.^2)).^((gamma+1)/(2*(gamma-1)))) - E;
M_e = fzero(fcn, j); % M_e = exit Mach
while M_e < 1
    M_e = fzero(fcn, j+1);
end
p_ratio=(1+(((gamma-1)/2)*(M_e^2)))^(gamma/(gamma-1)); % pressure ratio pc/pe

%plot perimeter lines
[x_slice, y_slice]=f_perimeter(p, r_rocket);

%star point slice starting nodes 
[x0, y0, r0, theta0, A0, pc0, rr0]=f_initial_slice(n_tip, n_fillet, n_side, p, eps, r_tip, r_fillet, z, h, A_t, rho_b, c_star, a, n);

%star point slice nodes for all iterations and relative chamber pressures
[x, y, r, theta, out_counter, yy]=f_slice_regression(x0, y0, r0, theta0, r_rocket, h, rr0, fps, p, A_t, A0, pc0, rho_b, c_star, a, n);

%cell arrays with star point slice coordinates (polar and cartesian) for all iterations
[cell_x, cell_y, cell_r, cell_theta, loops]=f_vectors_to_cells (x, y, r, out_counter);

%full grain coordinates for all iterations
[cell_x_total, cell_y_total, cell_r_total, cell_theta_total]=f_full_grain (loops, cell_r, cell_theta, p);

%re-calculates chamber pressures and burning areas for a smoother plot curve
[A_b, pc]=f_surfs_and_press_star (cell_x_total, cell_y_total, cell_r_total, r_rocket, h, A_t, rho_b, c_star, a, n);

%analytical calculations for burning areas during combustion
[A_b_an, yy_an, phase_B, phase_slivers, iteration_B, iteration_slivers]=f_analytical_method(loops, p, eps, th, w, f, l, h, A_t, rr0, pc0, a, n, rho_b, c_star, fps);
%------------------------------------------------------------------------%

%--------------------------GRAPHICS AND ANIMATIONS-----------------------%
%starting star point slice plot
figure(1)
hold on
plot([0 x_slice(end)],[0 y_slice(end)], 'k--');
plot(x_slice, y_slice, 'k', 'LineWidth', 2);
plot(x0, y0, 'r', 'LineWidth', 2);
xlim([0 r_rocket])
ylim([0 r_rocket])
title('Starting star point slice')
axis equal
% saveas(gcf,'Starting star point slice.png')
hold off

%star point slice regression plot
figure(2)
hold on
plot([0 x_slice(end)],[0 y_slice(end)], 'k--');
plot(x_slice, y_slice, 'k', 'LineWidth', 2);
plot(x0, y0, 'r');
xlim([0 r_rocket])
ylim([0 r_rocket])
daspect([1 1 1 ])
for i=2:loops
    plot(cell_x{i}, cell_y{i}, 'k')
end
title('Star point slice combustion surfaces')
axis equal
% saveas(gcf,'Star point slice combustion surfaces.png')
hold off

%full grain combustion surfaces plot
figure(3)
hold on
x_plot=cell_x_total{1};
y_plot=cell_y_total{1};
plot(x_plot,y_plot, 'r')
for i=2:loops
   x_plot=cell_x_total{i};
   y_plot=cell_y_total{i};
   if i<iteration_B
          plot(x_plot,y_plot, 'r')
   elseif i<iteration_slivers
          plot(x_plot,y_plot, 'b')
   else
          plot(x_plot,y_plot, 'm')
   end
end
theta_circ=linspace(-(pi/p) , (2*pi)-(pi/p) , size(x_plot,2));
r_circ=r_rocket*ones(1,size(x_plot,2));
x_circ=r_circ.*cos(theta_circ);
y_circ=r_circ.*sin(theta_circ);
plot(x_circ, y_circ, 'k', 'LineWidth', 2);
title('Full grain combustion surfaces')
axis equal
% saveas(gcf,'Full grain combustion surfaces.png')
hold off

%analytical vs numerical approach comparison plot
figure(4)
hold on
plot(yy, A_b, 'r', 'LineWidth', 2, 'DisplayName',' Numerical solution');
plot(yy_an, A_b_an, 'b', 'LineWidth', 1, 'DisplayName',' Analytical solution');
plot([phase_B phase_B],[0 1.1*max(A_b)], 'k--', 'DisplayName',' Phase B transition');
plot([phase_slivers phase_slivers],[0 1.1*max(A_b)], 'k', 'DisplayName',' Slivers transition');
xlim([0 1.1*yy(end)])
lgd=legend('Location', 'southwest');
lgd.FontSize=8;
title('Analytical vs numerical approach comparison');
xlabel('y [m]') 
ylabel('A_b [m^2]') 
grid on
% saveas(gcf,'Analytical vs numerical approach comparison.png')
hold off

%times vector initialisation
for i=1:loops
   time_vec(i)=i-1; 
end
time_vec=time_vec./fps;

%thrusts vector and toal impulse
CF=(Gam.*sqrt(((2*gamma)./(gamma-1)).*(1-((1./p_ratio).^((gamma-1)./gamma)))))+(E.*((1./p_ratio)-(p_amb./pc)));
F=A_t.*pc.*CF;
F(F<0)=0;
F(imag(F)~=0)=0;
I_t=trapz(time_vec,F);

%print on screen of some simulation results
disp(['Area ratio: ', num2str(E,4)]);
disp(['Exit Mach: ', num2str(M_e,4)]);
disp(['Characteristic velocity: ', num2str(c_star,6), ' m/s']);
disp(['Combustion duration: ', num2str(time_vec(end)), ' s']);
disp(['Total impulse: ', num2str(I_t,4), ' Ns']);

% %combustion animation frames
% [M_combustion]=f_combustion_animation (cell_x_total, cell_y_total, r_circ, theta_circ, loops, r_rocket);
% 
% %uncomment these lines if you want instant playback of the animation
% figure('visible','on');
% movie(M_combustion)
% 
% %saving video in current directory
% video = VideoWriter('combustion_video.avi');
% video.FrameRate = 10; % <-- to see in real time speed should set FrameRate=fps
% open(video)
% writeVideo(video,M_combustion)
% close(video)
% 
% %thrust-time curve animation frames
% [M_thrust]=f_thrust_time_animation (F, time_vec, loops);
% 
% %uncomment these lines if you want instant playback of the animation
% figure('visible','on');
% movie(M_thrust);
% 
% %saving video in current directory
% video = VideoWriter('thrust_time_curve_video.avi');
% video.FrameRate = 10; % <-- to see in real time speed should set FrameRate=fps
% open(video)
% writeVideo(video,M_thrust)
% close(video)
%------------------------------------------------------------------------%

toc