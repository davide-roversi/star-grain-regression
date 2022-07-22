function [cell_x_total, cell_y_total, cell_r_total, cell_theta_total]=f_full_grain (loops, cell_r, cell_theta, p)

delta_theta=(2*pi)/p; %angular increment for every star point slice of the full star

%IMPORTANT: "i" is the temporal iteration index
for i=1:loops
    r_tmp_up=cell_r{i}; %first slice upper nodes radial coordinates
    r_tmp_down=flip(r_tmp_up); %first slice bottom nodes radial coordinates
    theta_tmp_up=cell_theta{i}; %first slice upper nodes angular coordinates
    theta_tmp_down=-flip(theta_tmp_up); %first slice bottom nodes angular coordinates
    
    r_initial=[r_tmp_down, r_tmp_up]; %first slice all nodes radial coordinates
    theta_initial=[theta_tmp_down, theta_tmp_up]; %first slice all nodes angular coordinates
    
    r_final=r_initial; %initializing vector that will contain all radial coordinates for full grain at i-iteration
    theta_final=theta_initial; %initializing vector that will contain all angular coordinates for full grain at i-iteration
    
    %creation of radial and angular coordinates of full grain at i-iteration
    for j=1:p-1
       theta_tmp=theta_initial+(j*delta_theta);

       r_conc=[r_final r_initial];
       theta_conc=[theta_final theta_tmp];

       r_final=r_conc;
       theta_final=theta_conc;
    end
   
   %from polar coordinates to cartesian coordinates
   x_final=r_final.*cos(theta_final);
   y_final=r_final.*sin(theta_final);

   %filling i-row of cells that will contain coordinates of every node at every temporal iteration
   cell_r_total{i}=r_final;
   cell_theta_total{i}=theta_final;
   cell_x_total{i}=x_final;
   cell_y_total{i}=y_final;
end
end