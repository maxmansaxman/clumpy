%this is a script to read and process all the data on pressure baselines
%during a week, using an ETH-syle pressure baseline correction

%also the first crack at scripting the reading and importing of files
close all
clear all

cd ~/Dropbox/CarbonateClumping/PressureBaselineCorrection/Max_April_2014

files=dir('*.csv'); %taking only the csv data we want  
len=size(files,1);
N=cell(len,2);

for i=1:len
    fileName=files(i).name;
    A=importdata(fileName);
    fileNameTrim=fileName(1:end-4);
    N{i,1}=array2table(A.data);
    N{i,2}=fileNameTrim;
    N{i,1}.Properties.VariableNames={'num','HV','mass44','mass45',...
        'mass46','mass47','mass48','mass49'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


