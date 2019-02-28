close all;
path='/home/chen/Desktop/project/SCFT_results/O2/'
files=dir([path 'solution*.txt']);
result=cell(length(files),1);
detailedResult=cell(length(files),1);
for i=1:length(files)
    result{i}.name=files(i).name;
    fileID = fopen([path files(i).name],'r');
    tline = fgets(fileID);
    tmp = sscanf(tline, 'N= %d, ERROR= %e',2);
    result{i}.N=tmp(1);
    result{i}.error=tmp(2);
    tline = fgets(fileID);
    mean_field_free_energy = sscanf(tline, '%*s %f',1);
    result{i}.mean_field_free_energy=mean_field_free_energy;
    [A,cnt] = fscanf(fileID,'%*f,%f,%f',[2 result{i}.N]);
    A=A.';
    result{i}.x=A(:,1);
    result{i}.theta=A(:,2);
    fclose(fileID);
end

%% read Detailed Solution
files=dir([path 'detailed*.txt']);
for i=1:length(files)
    detailedResult{i}.name=files(i).name;
    fileID = fopen([path files(i).name],'r');
    detailedResult{i}.name=files(i).name;
    fileID = fopen([path files(i).name],'r');
    tline = fgets(fileID);
    tmp = sscanf(tline, 'N= %d, ERROR= %e',2);
    detailedResult{i}.N=tmp(1);
    detailedResult{i}.error=tmp(2);
    tline = fgets(fileID);
    mean_field_free_energy = sscanf(tline, '%*s %f',1);
    detailedResult{i}.mean_field_free_energy=mean_field_free_energy;
    [A,cnt] = fscanf(fileID,'%*f,%f,%f',[2 2^18]);
    A=A.';
    detailedResult{i}.x=A(:,1);
    detailedResult{i}.theta=A(:,2);
    fclose(fileID);
 end

%% sort files based on N
A=struct2cell(cell2mat(result));
B=cell2mat(A(2,:).');
[~, order] = sort(B);
result = result(order, :);
detailedResult = detailedResult(order, :);
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
xlabel('m');

%% plot detailedResult
figure;
for i=1:length(detailedResult)
    x=detailedResult{i}.x;
    theta=detailedResult{i}.theta;
    name=['Detailed N=' num2str(result{i}.N)];
    semilogx (x,theta,'*-','DisplayName',name);
    hold on;
    xlim([0 1]);
    legend show
end


