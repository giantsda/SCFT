function read_input(filename)
%% regular my input file
% global info;
% fprintf('reading %s...\n',filename);
% fid = fopen(filename,'r');
% if( fid==-1 )
%     error('Can''t open the file.');
% end
% str = fgets(fid);
% N = sscanf(str,'%*s %d', 1);
% fgets(fid); % skip x=0;
% while ~feof(fid)
%     str = fgets(fid);
%     data = sscanf(str,'%d,%f,%f', 3);
%     mesh(data(1))=data(2);
%     yita_middle_1D(data(1))=data(3);
% end
% fclose(fid);
% yita_middle_1D(end)=[];
% yita_middle_1D=yita_middle_1D.';
% mesh=[0 mesh];
% info.N=N;
% info.mesh=mesh;
% info.yita_middle_1D=yita_middle_1D;


%% Dr.Wang's results file
global info;
fprintf('reading %s...\n',filename);
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
end
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
N = sscanf(str,'m = %d', 1)+1;
for i=1:7
    fgets(fid); % skip head;
end
i=1;
while ~feof(fid)
    str = fgets(fid);
    data = sscanf(str,'%f %f %f', 3);
    mesh(i)=data(1)*info.L;
    yita_middle_1D(i)=data(3);
    i=i+1;
end
fclose(fid);
yita_middle_1D(end)=[];
yita_middle_1D=yita_middle_1D.';
mesh=[0 mesh];
info.N=N;
info.mesh=mesh;
info.yita_middle_1D=yita_middle_1D;
