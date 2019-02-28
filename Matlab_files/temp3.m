% close all;
% 
% l=3.723743573321600;
% x=Expm1024n2048IE(:,1)*l;
% y=Expm1024n2048IE(:,3);
% semilogx (x,y,'*-','DisplayName',name);
% hold on;
% 
% xx=result{6}.x;
% theta=result{6}.theta;
% name=['Detailed N=' num2str(result{i}.N)];
% semilogx (xx,theta,'-','DisplayName',name);
% hold on;
% xlim([0.05 1]);
% legend show
% 
% haha=theta-y;

%% load results.txt

folder_path='/home/chen/Desktop/project/SCFTTemp/DEALII_SCFT/solution_store.txt' ;
data=[];
delimiterIn = ','; %read txt file
headerlinesIn=0;
data_m = importdata(folder_path,delimiterIn,headerlinesIn);