function list_materials()

[~,materials]=refractive_index('FS',1);

string='Included materials (with call handle):\n';

for i=1:length(materials)
    string_mat=strcat(materials{i},'\n');
    string=strcat(string, string_mat);
end

fprintf(string)


end

