function [M_combustion]=f_combustion_animation (cell_x_total, cell_y_total, r_circ, theta_circ, loops, r_rocket)

x_circ=r_circ.*cos(theta_circ); %x coordinates of cylindrical combustion chamber perimeter nodes
y_circ=r_circ.*sin(theta_circ); %y coordinates of cylindrical combustion chamber perimeter nodes
x_plot=cell_x_total{1}; %x coordinates of first grain surface nodes (t=0)
y_plot=cell_y_total{1}; %y coordinates of first grain surface nodes (t=0)
M_combustion(loops) = struct('cdata',[],'colormap',[]); %frames matrix initialization

figure('visible','off'); %disabling view of every frame
hold on
plot(x_plot, y_plot, 'r', 'LineWidth', 1) %grain geometry plot
plot(x_circ, y_circ, 'k', 'LineWidth', 3) %rocket combustion chamber perimeter plot
xlim([-r_rocket r_rocket])
ylim([-r_rocket r_rocket])
daspect([1 1 1 ])
axis equal
x_color=[x_plot, flip(x_circ)]; %x coordinates of both grain and comb chamb nodes (to use "color fill")
y_color=[y_plot, flip(y_circ)]; %y coordinates of both grain and comb chamb nodes (to use "color fill")
c=[1*ones(1,size(x_plot,2)), 0*ones(1,size(x_circ,2))]; %color gradient intensities vector
fill(x_color, y_color, c, 'EdgeColor', 'none'); %color fill grain portion
colormap(autumn) %colormap style
M_combustion(1)=getframe(gcf); %first frame in frames matrix
hold off

%repeating above steps for every iteration
for i=2:loops
   x_plot=cell_x_total{i};
   y_plot=cell_y_total{i};
   figure('visible','off');
   hold on
   plot(x_plot, y_plot, 'r', 'LineWidth', 1)
   plot(x_circ, y_circ, 'k', 'LineWidth', 3)
   xlim([-r_rocket r_rocket])
   ylim([-r_rocket r_rocket])
   daspect([1 1 1 ])
   axis equal
   x_color=[x_plot, flip(x_circ)];
   y_color=[y_plot, flip(y_circ)];
   c=[(1)*ones(1,size(x_plot,2)), 0*ones(1,size(x_circ,2))];
   fill(x_color, y_color, c, 'EdgeColor', 'none');
   colormap(autumn)
   M_combustion(i)=getframe(gcf);
   hold off
end

end
