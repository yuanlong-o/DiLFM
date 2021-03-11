function result = glob(directory, pattern)
% glob, mimic the function in python
d = fullfile(directory, pattern);
files = dir(d);

result = cell(numel(files), 1);
for i = 1:numel(result)
    result{i} = fullfile(directory, files(i).name);
end
