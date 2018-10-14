 global info;
fileID = fopen('solution_matlab_N=33','w');
fprintf(fileID,'N=33\n');
x_old=[0;x_old;0];
for i=1:length(x_old)
    fprintf(fileID,'%d,%2.15f,%2.15f\n',i-1,info.mesh(i),x_old(i));
end
fclose(fileID);