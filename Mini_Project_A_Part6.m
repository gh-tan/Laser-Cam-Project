% Guohao Tan, MAE 150, Mini Project A Part 6
% DYnamic Analysis of Cam
close all; clear; clc

format shortG
%% P6 a)
w = 2.1; % constant angular velocity in rad/s

% initialize angle and each angle interval
theta0 = 0; % in degree, same as following
theta1 = 120;
theta2 = 200;
theta3 = 315;
theta4 = 360;
theta1_gap = theta1 - theta0; % in degree, degree between each interval
theta2_gap = theta2 - theta1;
theta3_gap = theta3 - theta2;
theta4_gap = theta4 - theta3;

% initialize height displacement, y
L1 = 10 - 0; % in mm, same as the following
L2 = 0;
L3 = abs(10 - 0);
L4 = 0;

% 0 < theta < 60, harmonic rise from 0mm to 10mm
theta1_ = linspace(theta0,theta1,120*100); % in deg, 1 deg per interval
beta1 = deg2rad(theta1_gap); % beta in radius, its a constant
y1 = L1/2*(1-cos(pi/beta1*deg2rad(theta1_)));
v1 = pi/beta1*(L1/2)*w*sin(pi/beta1*deg2rad(theta1_)); % angle need to be in rad, not deg!
% a1 = L1/2*((pi*w/beta1).^2)*cos(pi/beta1*deg2rad(theta1_)); % this is just a point
% jerk1 = -L1/2*((pi*w/beta1)^3)*sin(pi*theta1_/beta1);


% 120 < theta < 200, dwell
theta2_ = linspace(theta1,theta2,(200-120)*100); % in deg
y2(1:length(theta2_)) = 10; % 10 mm
v2(1:length(theta2_)) = 0;
% a2(1:length(theta2_)) = 0;
% jerk2(1:length(theta2_)) = 0;

% 200 < theta < 315, linear fall
theta3_ = linspace(theta2-theta2,theta3-theta2,(315-200)*100); % in deg, with phase shift?
beta3 = deg2rad(theta3_gap); % in radius
y3 = L3-L3/beta3*deg2rad(theta3_);
v3(1:length(theta3_)) = L3/beta3*w;
% a3 = 0;
% jerk3 = 0;

% 315 < theta < 360, dwell
theta4_ = linspace(theta3,theta4,(360-315)*100); % in deg
y4(1:length(theta4_)) = 0; % 0 mm
v4(1:length(theta4_)) = 0;
% a4(1:length(theta4_)) = 0;
% jerk4(1:length(theta4_)) = 0;

% plot displacement vs angle
figure(1)
plot(theta1_,y1,'r')
axis([0 360,0 20])
hold on
grid on
title('Part 6) Displacement Vs. Theta for 1 rev')
xlabel('Angle, theta in degree');ylabel('Displacement, y in mm')
ylim([-2 12])
plot(theta2_,y2,'b')
plot(200+theta3_,y3,'g')
plot(theta4_,y4,'m')

% sum all theta and (velocity, displacement) of y
theta = deg2rad([theta1_ theta2_ 200+theta3_ theta4_]);
y = [y1 y2 y3 y4];
v = [v1 v2 v3 v4];



%% part b)

m = 0.5; % in kg
k = 150; % in N/m
DR = 0.5; % damping ratio

% DT = theta(end)/w * ones(length(theta),1);% constant DT, DT = Delta T, a matrix
DT = (2*pi/(360*100)) * ones(360*100,1);

dydt = v; % velocity for cam

% ICs
x = zeros(1,length(DT)); % make empty matrix for x, x(1) = 0 >> IC
dxdt = zeros(1,length(DT)); % dxdt(1) = 0 >> IC
dx2dt2 = zeros(1,length(DT));

% use for loop to loop through eqns for euler method
for ii = 1:36000-1 % Used 36000 points(theta) in part a)
    
    dx2dt2(ii) = 2.*DR./m.*sqrt(k*m).*dydt(ii)+k/m.*y(ii) ...
        -2*DR/m*sqrt(k*m).*dxdt(ii)-k/m*x(ii); 
    
    dxdt(ii+1) = dxdt(ii) + dx2dt2(ii).*DT(ii);
    x(ii+1) = x(ii) + dxdt(ii).*DT(ii);
    
    
end

plot(rad2deg(theta),x,'k')
legend('Hamonic rise','dwell','linear fall','dwell','BY Euler Method','location','best')


