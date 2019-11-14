% Data citation:
%  UCAR/NCAR - Earth Observing Laboratory. 2013. PCAPS ISFS 1 second data.
%  Version 1.0. UCAR/NCAR - Earth Observing Laboratory.
%  https://doi.org/10.5065/D6QV3JRP. Accessed 24 Jan 2018.

% Data is from Salt Lake basin, on 2010-11-10


close all
clear all
% Import and plot data original 1
data = csvread('wind_data3.csv');
figure
plot(data)
xlabel("time (s)")
ylabel("wind speed (m/s)")


% Compute mean by adding the speed at each time and dividing by the total
% number of measurements in the data set
sum=0;
for i=1:length(data)
    sum=sum+data(i);
end
average=sum/length(data)

% Variance=average(u'(t)^2)
% Compute variance by keeping running total of difference between the speed
% at a given time and the average, squared
varsum=0;
for i=1:length(data)
    varsum=varsum+(data(i)-average)^2;
end
variance=varsum/length(data) % Average over all times to find variance (assuming population)
st_dev=sqrt(variance) % Standard deviation is square root of variance

% Skewness=average(u'(t)^3)/sigma^3
% Compute skewness by finding the value of the numerator of the formula at
% each given time, summing, and dividing by the total number of
% measurements and by the denominator
skewsum=0;
for i=1:length(data)
    skewsum=skewsum+(data(i)-average)^3;
end
skew=skewsum/length(data)/st_dev^3

% Flatness=average(u'(t)^4)/sigma^4
% Compute flatness by finding the value of the numerator of the formula at
% each given time, summing, and dividing by the total number of
% measurements and by the denominator
flatsum=0;
for i=1:length(data)
    flatsum=flatsum+(data(i)-average)^4;
end
flatness=flatsum/length(data)/variance^2

% Create CDF by dividing the data into 1000 bins covering an equal range,
% then using the results to find the proportion of the data points that are
% below a given value.  This was done by finding the portion of the data
% points that were in each bin and adding this to the portion of data
% points that were lower in value than the minimum value of the bin. The
% speed corresponding to each value of the CDF is found by using the
% maximum value of the bin.
% The CDF was used instead of the PDF because it allowed for smaller bin
% sizes without being strongly affected by random variability.
[N,edges]=histcounts(data,1000); % divide data into bins
% Since the value of the CDF at a given speed depends on its value at a
% lower speed, the first term must be found separately.
F(1)=N(1)/length(data);
u(1)=edges(2);
for i=2:length(N) % Use a loop for the remainder of the terms
    F(i)=F(i-1)+N(i)/length(data);
    u(i)=edges(i+1);
end
figure
plot(u,F)
hold on

% Fit using Weibull distribution
% Weibull CDF: F(u)=1-exp(-(u/lambda)^K)
% Weibull PDF: f(u)=K(u^(k-1)/lambda^k)exp(-(u/lambda)^K)
% Write each of these as a function and combine the parameters in the Weib 
% array, where Weib(1)=lambda and Weib(2)=K
Func=@(Weib,xdata)(1-exp(-(xdata/Weib(1)).^Weib(2)));
func=@(Weib,xdata)Weib(2)*(xdata.^(Weib(2)-1)/Weib(1).^Weib(2)).*exp(-(xdata/Weib(1)).^Weib(2));
Weib0=[1.3 2]; % Give starting values of the parameters
% Fit the CDF of the Weibull distribution to the data using a least squares
% curve fit
[Weib,resnorm,~,exitflag,output] = lsqcurvefit(Func,Weib0,u,F);
Weib
plot(u,Func(Weib,u))

% Plot histogram with Weibull PDF
% Make the two scales match by making the maximum y-value of the PDF 1 and
% the maximum y-value of the histogram equal the number of measurements
% multiplied by the width of a bin.
figure
yyaxis left
histogram(data)
ylim([0,0.1*length(data)]);
ylabel("frequency")
xlabel("speed (m/s)")
title("Wind Speed Distribution")
yyaxis right
plot(u,func(Weib,u))
ylim([0,1]);
ylabel("probability")

% Find the autocorrelation for each time delay from 1 second to 1850 seconds
% (where it crosses zero).  Beyond 1850 seconds, the values become too 
% erratic to effectively use.
% Uses the function at the bottom of this script
for i=1:1850
    auto(i)=autocorrelation(i,data,average); % formula from class notes
    autodl(i)=auto(i)/variance; % dimensionless version
end
figure
plot(autodl)
grid on
figure
semilogx(autodl)
xlabel("tau")
ylabel("Ruu")
grid on
title("Time Delau Autocorrelation")

% Compute integral time scale by adding the dimensionless value of the
% autocorrelation at each tau.
its=0;
for i=1:length(autodl)
    its=its+autodl(i);
end
its

% Find the energy spectrum by taking the Fourier transform of the
% autocorrelation data
figure
Y=abs(fft(auto)/length(auto));
f=(0:(length(auto))-1)/length(auto);
%  Y2=Y(1:length(auto)/2+1);
%  f2 = (0:(length(auto)/2))/length(auto);
%  plot(f2(2:end),Y2(2:end))
%  figure
plot(f,Y)
ylabel("Energy (m^2/s)")
xlabel("Frequency (Hz)")
title("Kinetic Energy Spectrum")
grid on



% Autocorrelation function for given tau
function Ruu=autocorrelation(tau,data,average)
i=1; % variable to iterate
n=0; % Keep track of how many numbers are being averaged
ruusum=0; 
% Use this rather than the standard for loop to ensure that calculations
% only occur for times in which u'(t) and u'(t+tau) are both available
while (i+tau)<=length(data)
    n=n+1;
    % add value at each time to sum
    ruusum=ruusum+(data(i+tau)-average)*(data(i)-average); 
    i=i+1;
end
Ruu=ruusum/n;
end