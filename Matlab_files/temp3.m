folder_path='C:\Users\chenshen\Downloads' ;
file_name='Exp_m32_n2048_IE.res';
data=[];
filename1= [folder_path '\' file_name ] ;
disp(filename1);
delimiterIn = ' '; %read txt file
headerlinesIn = 9;
data_m = importdata(filename1,delimiterIn,headerlinesIn);
data=data_m.data;
yita_answer=data(:,3);
