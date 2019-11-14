close all
clear all
filename = 'wind data qwer.csv';
num = csvread('2010_MET_DATA_noHeaders.csv');
speedranges=[0 5 10 15 20 25];
[figure_handle,count,speeds,directions,Table] = WindRose(num(:,7),num(:,6),'vWinds',speedranges);
[figure_handle,count,speeds,directions,Table] = WindRose(num(:,12),num(:,11));
[figure_handle,count,speeds,directions,Table] = WindRose(num(:,17),num(:,16));
[figure_handle,count,speeds,directions,Table] = WindRose(num(:,22),num(:,21));
[figure_handle,count,speeds,directions,Table] = WindRose(num(:,27),num(:,26));
figure
histogram(num(:,7))