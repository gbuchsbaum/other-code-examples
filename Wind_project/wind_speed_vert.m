close all
clear all

% Import wind speed and direction data
filename = 'wind data qwer.csv';
num = csvread('data_2015.csv');

% Limit to data in prevailing wind direction.  Frome wind rose, directions
% between 180 degrees and 210 degrees were used.
% For each time instance in which the wind is coming from the correct
% direction at 80 m, save the speed and the speed cubed at each height.
j=1;
for i=1:length(num)
    if num(i,7)>180
        if num(i,7)<210
            s80(j,1)=num(i,6);
            s80(j,2)=num(i,6).^3;
            s90(j,1)=num(i,11);
            s90(j,2)=num(i,11).^3;
            s100(j,1)=num(i,16);
            s100(j,2)=num(i,16).^3;
            s110(j,1)=num(i,21);
            s110(j,2)=num(i,21).^3;
            s120(j,1)=num(i,26);
            s120(j,2)=num(i,26).^3;
            j=j+1;
        end
    end
end

% Find the average wind speed at each height and save it in an array that
% includes the height.  Also find the average of the cubed wind speeds and
% take the cube root of that value to represent the average power
% available.
avg(:,1)=80:10:120;
avgpow(:,1)=80:10:120;
avg(1,2)=mean(s80(:,1));
avgpow(1,2)=(mean(s80(:,2)))^(1/3);
avg(2,2)=mean(s90(:,1));
avgpow(2,2)=(mean(s90(:,2)))^(1/3);
avg(3,2)=mean(s100(:,1));
avgpow(3,2)=(mean(s100(:,2)))^(1/3);
avg(4,2)=mean(s110(:,1));
avgpow(4,2)=(mean(s110(:,2)))^(1/3);
avg(5,2)=mean(s120(:,1));
avgpow(5,2)=(mean(s120(:,2)))^(1/3);

% Plot the average wind speed at each height using both methods.
scatter(avg(:,1),avg(:,2))
figure
scatter(avgpow(:,1),avgpow(:,2))

% Create fit.  Make a vector with just the height, a vector with just
% the average speeds, and a vector with just the speed corresponding to the
% average power availability.
z=(avg(:,1));
u=(avg(:,2));
up=(avgpow(:,2));

% Create a fit type using the power law vertical profile from class notes.
% This uses the height as the independent variable and the speed as the
% dependent variable, and has coefficients representing the reference
% height zr, the speed at the reference height uzr, and the Hellmann
% exponent a.  Initial values and lower limits are given to add ease to
% help the fit converge.
myfittype = fittype('uzr*(z/zr)^a',...
    'dependent',{'u'},'independent',{'z'},...
    'coefficients',{'a','uzr','zr'});
options = fitoptions(myfittype);
options.StartPoint = [0.2 10 6];
options.Lower = [0 0 0];

myfittype2 = fittype('ust/K*log(z/z0)',...
    'dependent',{'u'},'independent',{'z'},...
    'coefficients',{'K','ust','z0'});
options = fitoptions(myfittype);
options.StartPoint = [0.1 10 6];
options.Lower = [0 0 0];

% Calculate the fits for both methods of determining average speed, and
% plot them with the data to confirm that the fit is accurate.
myfit = fit(z,u,myfittype,options)
myfit2 = fit(z,up,myfittype2,options)
figure
plot(myfit,z,u)
myfitpow = fit(z,up,myfittype,options)
figure
plot(myfitpow,z,up)
figure
plot(myfit2,z,up)

% Create an array of the estimated wind speed at heights from 50m to 150m,
% using the fit identified.  This only uses the power-based average.
% Finally, save this array to a .mat file so it can be used as part of the
% optimization.
speedfit(:,1)=0:300;
speedfit(:,2)=myfitpow(speedfit(:,1));
figure
plot(speedfit(:,1),speedfit(:,2))
speedfit2(:,1)=0:300;
speedfit2(:,2)=myfit2(speedfit2(:,1));
figure
plot(speedfit(:,1),speedfit(:,2))
hold on
plot(speedfit2(:,1),speedfit2(:,2))
legend("log","power")
save('speedfit.mat','speedfit')