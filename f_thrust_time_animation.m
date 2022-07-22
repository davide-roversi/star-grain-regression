function [M_thrust]=f_thrust_time_animation (F, time_vec, loops)

M_thrust(loops) = struct('cdata',[],'colormap',[]); %frames matrix initialization

figure('visible','off'); %disabling view of every frame
plot(time_vec(1,1), F(1,1), 'r', 'LineWidth', 2); %first point of curve
xlim([0 1.1*time_vec(1,end)])
ylim([0 1.1*max(F)])
title('Thrust-time curve');
xlabel('t [s]')
ylabel('Thrust [N]')
grid on
M_thrust(1)=getframe(gcf); %first frame in frames matrix

%repeating above steps for every iteration
for i=2:loops
    figure('visible','off');
    plot(time_vec(1,[1:i]), F(1,[1:i]), 'r', 'LineWidth', 2);
    xlim([0 1.1*time_vec(1,end)]);
    ylim([0 1.1*max(F)]);
    title('Thrust-time curve');
    xlabel('t [s]')
    ylabel('Thrust [N]')
    grid on
    M_thrust(i)=getframe(gcf);
%     if i==loops
%        saveas(gcf,'Thrust-time curve.png') 
%     end
end

end
