%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cell_data=matrix2cell(matrix)

dim_matrix=size(matrix);
cell_data=cell(dim_matrix(1),1);


for i=1:dim_matrix(1)
    temp_matrix=squeeze(matrix(i,:,:));
    cell_data{i}=temp_matrix;
end


end
