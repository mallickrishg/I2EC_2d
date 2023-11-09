function filename = filenameGenerator(foldername)
% function to read a results folder and create a new sub-folder to store
% results.mat
% Rishav Mallick, JPL 2023

fout = dir([foldername '/0*']);

num = length(fout)+1;

if num<10
    outputfolder = [foldername '/00000' num2str(num)];
elseif num>=10 & num<100
    outputfolder = [foldername '/0000' num2str(num)];
elseif num>=100 & num<1000
    outputfolder = [foldername '/000' num2str(num)];
elseif num>=1000 & num<10000
    outputfolder = [foldername '/00' num2str(num)];
else
    error('too many experiments, create a new folder!')
end

mkdir(outputfolder)
filename = [outputfolder '/results.mat'];

end