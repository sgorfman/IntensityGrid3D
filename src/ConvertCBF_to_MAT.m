function ConvertCBF_to_MAT
%It will go through all the cbf files in the folder and copies all of them to matlab format
d = dir('*.cbf'); %Find all the cbf files 
N = numel(d); 
for n = 1:N
    FileName_n = d(n).name
    F = cbfread(FileName_n);
    [~,FileName_n_new] = fileparts(FileName_n);
    FileName_n_new = [FileName_n_new,'.mat'];
    save(FileName_n_new,'F');
    delete(FileName_n);
end