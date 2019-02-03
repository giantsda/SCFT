close all;
path='../DEALII_SCFT/'
files=dir([path 'solution*.txt']);
result=cell(length(files),1);
for i=1:length(files)
    result{i}.name=files(i).name;
    fileID = fopen([path files(i).name],'r');
    tline = fgets(fileID);
    N = sscanf(tline, '%*s %d',1);
    result{i}.N=N;
    [A,cnt] = fscanf(fileID,'%*f,%f,%f',[2 N]);
    A=A.';
    result{i}.x=A(:,1);
    result{i}.theta=A(:,2);
    tline = fgets(fileID);
    tline = fgets(fileID);
    mean_field_free_energy = sscanf(tline, '%*s %f',1);
    result{i}.mean_field_free_energy=mean_field_free_energy;
    fclose(fileID);
end

%% sort files based on N
A=struct2cell(cell2mat(result));
B=cell2mat(A(2,:).');
[~, order] = sort(B);
result = result(order, :)
A=struct2cell(cell2mat(result));
B=cell2mat(A(2,:).');

mean_field_free_energy_s=[];
for i=1:length(result)
    x=result{i}.x;
    theta=result{i}.theta;
    name=['N=' num2str(result{i}.N)];
    semilogx (x,theta,'*-','DisplayName',name);
    hold on;
    xlim([0 1]);
    legend show
    mean_field_free_energy_s=[mean_field_free_energy_s result{i}.mean_field_free_energy];
end

figure;
absDiff=abs(mean_field_free_energy_s-mean_field_free_energy_s(end));
loglog (B,absDiff,'*-');
title('|Fm-F1024|');
xlabel('m')
 