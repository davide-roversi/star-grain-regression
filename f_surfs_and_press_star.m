function [A_b, pc]=f_surfs_and_press_star (cell_x_total, cell_y_total, cell_r_total, r_rocket, h, A_t, rho_b, c_star, a, n)

for i=1:size(cell_x_total,2)
    S_b=0; %perimeter initialization
    x_count=cell_x_total{i}; %extracting only the x coordinates to be considered
    y_count=cell_y_total{i}; %extracting only the y coordinates to be considered
    r_count=cell_r_total{i}; %extracting only the radial coordinates to be considered
    for j=2:size(x_count,2) %calculating combustion perimeter at i-iteration
        if r_count(1,j)< r_rocket
           S_b=S_b+sqrt(((x_count(j)-x_count(j-1))^2)+((y_count(j)-y_count(j-1))^2));
        end
    end
    A_b(i)=S_b*h; %combustion surface at i-iteration
    K=A_b(i)/A_t; %Klemmung factor at i-iteration
    pc(i)=(K*rho_b*c_star*a)^(1/(1-n)); % [Pa] chamber pressure at i-iteration
end

end