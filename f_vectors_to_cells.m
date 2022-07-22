function [cell_x, cell_y, cell_r, cell_theta, loops]=f_vectors_to_cells (x, y, r, out_counter)

%creating cell arrays with polar and cartesian coordinates of every node of every star point slice surface in time

loops = size(r,1); %number of total iterations from combustion start to finish
cell_x{1}=x(1,:); 
cell_y{1}=y(1,:);
x_tmp=cell_x{1};
y_tmp=cell_y{1};
cell_r{1}=sqrt((x_tmp.^2)+(y_tmp.^2)); 
cell_theta{1}=atan(y_tmp./x_tmp);
for i=2:loops
    in_vec=find(out_counter(i,:)==0);
    cell_x{i}=x(i,in_vec); 
    cell_y{i}=y(i,in_vec);
    x_tmp=cell_x{i};
    y_tmp=cell_y{i};
    cell_r{i}=sqrt((x_tmp.^2)+(y_tmp.^2)); 
    cell_theta{i}=atan(y_tmp./x_tmp);
end

end
