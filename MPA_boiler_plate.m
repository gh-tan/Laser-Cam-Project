%% Initialize
clear;
clc;
close all;

%% Define some constants

r_f = 0.275;    % inches radius of follower
L = 2.5;     % inches cam separation difference from follower
R_b = 2.0;      % inches radius of base circle

%% Hand trace your pattern

figure(1);
clf;
title('Trace the desired cam profile', 'fontsize', 18);
axis([-1, 1, -1, 1]);

xx = zeros(100,1);
yy = zeros(100,1);

[x_left, y_left] = pol2cart(linspace(0,2*pi,100)', R_b*ones(100,1));
plot(x_left, y_left);
hold on;
plot(0,0,'ro')
plot(x_left + L, y_left);
plot(L,0,'bo')
axis equal;
axis([0, 2.5, 0.5, 2.5]);

plot([1, 1, 1.55, 1.55, 1], [1.55, 2.10, 2.10, 1.55, 1.55], '--k');

for ii=1:100
    
    title(['Trace the desired cam profile in the dashed box: ' num2str(ii)], 'fontsize', 18);
    
    [x_select, y_select, inp] = ginput(1);
    
    if ~isempty(inp)
       xx(ii) = x_select;
       yy(ii) = y_select;
        
    else
        break;
    end
    clf;
    plot(x_left, y_left);
    hold on;
    plot(0,0,'ro')
    plot(x_left + L, y_left);
    plot(L,0,'bo')
    plot([1, 1, 1.55, 1.55, 1], [1.55, 2.10, 2.10, 1.55, 1.55], '--k');
    axis equal;
    axis([0, 2.5, 0.5, 2.5]);
    
    plot(xx(1:ii), yy(1:ii), 'o-k');
end

xx = xx(1:ii-1);
yy = yy(1:ii-1);
xx(end+1) = xx(1);    % make sure is periodic
yy(end+1) = yy(1);    % make sure is periodic


plot(xx, yy, 'o-k');

%% save list of x,y points into text file
dlmwrite('x_y_data.txt',[xx(1:ii) yy(1:ii)],',');

%% plot x_i,y_i versus i
data_xy = dlmread('x_y_data.txt',',');
x_i = data_xy(:,1);
y_i = data_xy(:,2);
num_i = linspace(1,length(x_i),length(x_i));
figure(2)
plot(num_i,x_i,'ob')
title('plot of x_{i} and y_{i} Vs. i')
xlabel('i');ylabel('position point')
grid on
hold on
plot(num_i,y_i,'or')
legend('x_{i}','y_{i}','location','best')

% % draw the shape of alphabet
% figure(3)
% plot(x_i,y_i,'o')
% title('shape of assigned alphabet, Lambda(upper letter)')
% xlabel('x');ylabel('y')
% grid on
% axis equal
% xlim([1 1.5])

% here x_i = xx;y_i = yy so that we can just read date from files directly
% w/o redraw
xx = x_i;
yy = y_i;

%% Solve for alpha_i, beta_i, A_i, B_i, 


%%%%%%%%%%% Homework assignment

% write out eqn for alpha_i,beta_i,A_i,B_i
alpha_i = atan(yy./xx);
beta_i = atan(yy./(L-xx));
A_i = (xx./cos(alpha_i)) - r_f;
B_i = (yy./sin(beta_i)) - r_f;


%% Generata dummy wt


%%%%%%%%%%% Homework assignment
wt_i = linspace(0,2*pi,length(xx)); % range from 0 to 360

% for animation
wt = wt_i';

%% Solve for theta_left and theta_right


%%%%%%%%%%% Homework assignment
theta_i_L = (2*pi - wt_i') + alpha_i;
theta_i_R = (wt_i' + pi) - beta_i;


%% Generate the xy points of the cams


%%%%%%%%%%% Homework assignment
% write a matlab func to return the polar coord of left&right cam profile
% from (xi,yi) coord 

% polarFUN_Ai = @(xx,yy) A_i;
% % A_i_polar = polarFUN_Ai(xx,yy);
% 
% polarFUN_angle_L = @(xx,yy) theta_i_L;
% % theta_i_L_polar = polarFUN_angle_L(xx,yy);
% 
% polarFUN_Bi = @(xx,yy) B_i;
% % B_i_polar = polarFUN_Bi(xx,yy);
% 
% polarFUN_angle_R = @(xx,yy) theta_i_R;
% % theta_i_R_polar = polarFUN_angle_R(xx,yy);

% if change polar to cart coord
% [x_cam_left, y_cam_left] = pol2cart(theta_i_L,A_i);
% [x_cam_right, y_cam_right] = pol2cart(theta_i_R,B_i);

x_cam_left = A_i.*cos(theta_i_L);
y_cam_left = A_i.*sin(theta_i_L);
x_cam_right = B_i.*cos(theta_i_R);
y_cam_right = B_i.*sin(theta_i_R);


%% Draw the cam shapes 


%%%%%%%%%%% Homework assignment
% figure(3)
% polarplot(theta_i_L_polar,A_i_polar)
% figure(4)
% polarplot(theta_i_L_polar,B_i_polar)
figure(3)
plot(x_cam_left, y_cam_left,'b','linewidth',2)
title('Part2: 5) Cam Profiles')
xlabel('x coordinate');ylabel('y coordinate')
axis equal
grid on
hold on
plot(L+x_cam_right, y_cam_right,'r','linewidth',2) % need to add up the distance L for proper displayment of right cam
legend('Left Cam','Right Cam','location','best')

%% Part2 6) calculate the radius of curvature

% for left cam
% Dheta_theta_L = (theta_i_L(ii+1) - theta_i_L(ii));

% Ss_L = radius of curvature for left cam
% Ss_R = radius of curvature for right cam

% for ii = 2:(length(xx)-1)
%     ri_prime_L = (A_i(ii+1) - A_i(ii))./(theta_i_L(ii+1) - theta_i_L(ii));
%     ri_2_prime_L = (A_i(ii+1) -2*A_i(ii) + A_i(ii-1))./ ((theta_i_L(ii+1) - theta_i_L(ii)).^2);
%     Ss_L(ii) = (( (A_i(ii).^2) + (ri_prime_L.^2) ).^(3/2))./((A_i(ii).^2)+2.*(ri_prime_L.^2)-A_i(ii).*ri_2_prime_L);
%     if ii == length(xx)-1 % radius of curvature from cam end point to start point
%         ri_prime_L = (A_i(ii+1) - A_i(1))./(theta_i_L(ii+1) - theta_i_L(1));
%         ri_2_prime_L = (A_i(ii+1) -2*A_i(1) + A_i(ii-1))./ ((theta_i_L(ii+1) - theta_i_L(1)).^2);
%         Ss_L(length(xx)) = (( (A_i(1).^2) + (ri_prime_L.^2) ).^(3/2))./((A_i(1).^2)+2.*(ri_prime_L.^2)-A_i(1).*ri_2_prime_L);
%     end
% end

% 1st point to 2nd point
ri_prime_L = (A_i(2)+r_f - (A_i(1)+r_f))./(theta_i_L(2) - theta_i_L(1));
ri_2_prime_L = (A_i(2)+r_f -2*(A_i(1)+r_f) + (A_i(end)+r_f))./ ((theta_i_L(2) - theta_i_L(1)).^2);
Ss_L(1) = (( ((A_i(1)+r_f).^2) + (ri_prime_L.^2) ).^(3/2))./(((A_i(1)+r_f).^2)+2.*(ri_prime_L.^2)-(A_i(1)+r_f).*ri_2_prime_L);

for ii = 2:(length(xx)-1)
    ri_prime_L = (A_i(ii+1)+r_f - (A_i(ii)+r_f))./(theta_i_L(ii+1) - theta_i_L(ii));
    ri_2_prime_L = (A_i(ii+1)+r_f -2*(A_i(ii)+r_f) + (A_i(ii-1)+r_f))./ ((theta_i_L(ii+1) - theta_i_L(ii)).^2);
    Ss_L(ii) = (( ((A_i(ii)+r_f).^2) + (ri_prime_L.^2) ).^(3/2))./(((A_i(ii)+r_f).^2)+2.*(ri_prime_L.^2)-(A_i(ii)+r_f).*ri_2_prime_L);
    if ii == length(xx)-1 % radius of curvature from cam end point to start point, ii+1 = end point of ri (which is A_i+r_f)
        ri_prime_L = (A_i(1)+r_f - (A_i(ii+1)+r_f))./(theta_i_L(1) - theta_i_L(ii+1));
        ri_2_prime_L = (A_i(1)+r_f -2*(A_i(ii+1)+r_f) + (A_i(ii)+r_f))./ ((theta_i_L(1) - theta_i_L(ii+1)).^2);
        Ss_L(ii+1) = (( ((A_i(ii+1)+r_f).^2) + (ri_prime_L.^2) ).^(3/2))./(((A_i(ii+1)+r_f).^2)+2.*(ri_prime_L.^2)-(A_i(ii+1)+r_f).*ri_2_prime_L);
    end
end

figure(4)
theta_x_axis = linspace(0,2*pi,length(xx));
plot(theta_x_axis,Ss_L,'b')
title('R_{curvatureL} vs. Theta')
xlabel('Theta in radian');ylabel('R_{curvatureL} in inches');
grid on

% for right cam
% Dheta_theta_R = (theta_i_R(ii+1) - theta_i_R(ii));

% 1st point to 2nd point
ri_prime_R = (B_i(2)+r_f - (B_i(1)+r_f))./(theta_i_R(2) - theta_i_R(1));
ri_2_prime_R = (B_i(2)+r_f -2*(B_i(1)+r_f) + (B_i(end)+r_f))./ ((theta_i_R(2) - theta_i_R(1)).^2);
Ss_R(1) = (( ((B_i(1)+r_f).^2) + (ri_prime_R.^2) ).^(3/2))./(((B_i(1)+r_f).^2)+2.*(ri_prime_R.^2)-(B_i(1)+r_f).*ri_2_prime_R);

for ii = 2:(length(xx)-1)
    ri_prime_R = (B_i(ii+1)+r_f - (B_i(ii)+r_f))./(theta_i_R(ii+1) - theta_i_R(ii));
    ri_2_prime_R = (B_i(ii+1)+r_f -2*(B_i(ii)+r_f) + (B_i(ii-1)+r_f))./ ((theta_i_R(ii+1) - theta_i_R(ii)).^2);
    Ss_R(ii) = ((((B_i(ii)+r_f).^2)+(ri_prime_R).^2).^(3/2))./(((B_i(ii)+r_f).^2)+2.*((ri_prime_R).^2)-(B_i(ii)+r_f).*ri_2_prime_R);
    if ii == length(xx)-1 % radius of curvature from cam end point to start point
        ri_prime_R = (B_i(1)+r_f - (B_i(ii+1)+r_f))./(theta_i_R(1) - theta_i_R(ii+1));
        ri_2_prime_R = (B_i(1)+r_f -2*(B_i(ii+1)+r_f) + (B_i(ii)+r_f))./ ((theta_i_R(1) - theta_i_R(ii+1)).^2);
        Ss_R(ii+1) = ((((B_i(ii+1)+r_f).^2)+(ri_prime_R).^2).^(3/2))./(((B_i(ii+1)+r_f).^2)+2.*((ri_prime_R).^2)-(B_i(ii+1)+r_f).*ri_2_prime_R);
    end
    
end

figure(5)
theta_x_axis = linspace(0,2*pi,length(xx));
plot(theta_x_axis,Ss_R,'r')
title('R_{curvatureR} vs. Theta')
xlabel('Theta in radian');ylabel('R_{curvatureR} in inches');
grid on

figure(6)
theta_x_axis_L = linspace(0,2*pi,length(xx));
plot(theta_x_axis_L,Ss_L,'b')
hold on
title('R_{curvature L&R} vs. Theta')
xlabel('Theta in radian');ylabel('R_{curvature} in inches');
grid on
theta_x_axis_R = linspace(0,2*pi,length(xx));
plot(theta_x_axis_R,Ss_R,'r')
legend('R_{curvature L}','R_{curvature R}','location','best')

%% part2 7)
disp('')

%% part2 9 calculate the dynamic factor
m_laser = 0.102; % in kg
k_spring = 17.38; % in N/m
RPM = 90; % rev/mins
w_RPM = RPM*2*pi/60; % angular velocity, radians/second

dynamic_factor = m_laser/k_spring*(w_RPM^2);
fprintf('\nThe dynamic factor is %d \n',dynamic_factor)

%% make the follower and cams into homogeneous coordinates
% Note: This requires you have x_cam_left, y_cam_left, x_cam_right,
% y_cam_right vectors

[x_follower, y_follower] = pol2cart(linspace(0,2*pi,100)', r_f*ones(100,1));

follower(1,:) = x_follower;
follower(2,:) = y_follower;
follower(3,:) = 1;


% This is an arbitrary (2D) rotation matrix
% redefined for clockwise motion
R = @(theta) [cos(-theta)  sin(-theta) 0;
              -sin(-theta) cos(-theta) 0;
              0           0          1];

cam_left = [];
cam_left(1,:) = x_cam_left;
cam_left(2,:) = y_cam_left;
cam_left(3,:) = 1;

cam_left_coords(1,1:2) = [0, .5];
cam_left_coords(2,1:2) = [0, 0];
cam_left_coords(3,1:2) = [1, 1];

cam_right = [];
cam_right(1,:) = x_cam_right;
cam_right(2,:) = y_cam_right;
cam_right(3,:) = 1;

cam_right_coords(1,1:2) = [0, .5];
cam_right_coords(2,1:2) = [0, 0];
cam_right_coords(3,1:2) = [1, 1];


%% Now make an animation and plot the y, vel, and accel

figure(7)
set(gcf, 'color', 'w');

for i=1:length(A_i);
    
    ii = mod(i-1, length(A_i))+1;
    
    new_cam_left = R(wt(ii))*cam_left;
    new_cam_right = R(-wt(ii))*cam_right;
    
    new_coord_left = R(wt(ii))*cam_left_coords;
    new_coord_right = R(-wt(ii))*cam_right_coords;
    
    clf;
    plot(new_cam_left(1,:), new_cam_left(2,:), 'linewidth', 3);
    hold on;
    plot(new_coord_left(1,:), new_coord_left(2,:), 'k', 'linewidth', 3);
%     plot(x_left, y_left);
    plot(0,0,'ro')
    
    plot(new_cam_right(1,:) + L, new_cam_right(2,:), 'linewidth', 3);
    plot(new_coord_right(1,:) + L, new_coord_right(2,:), 'k','linewidth', 3);
%     plot(x_left + L, y_left);
    plot(L,0,'bo')
    
    axis equal;
    axis([0, 2.5, 0.5, 2.5]);
    
    plot(follower(1,:)+xx(ii), follower(2,:)+yy(ii), 'linewidth', 3);
    plot(xx(1:ii), yy(1:ii), '-k');
    plot(xx(ii), yy(ii), 'or');

    
    drawnow;
    
    % Uncomment these lines if you want to pause motion or save images
%     pause;
%     print([num2str(1000 + i) '.png'], '-dpng')
end

% If you have imagemagick installed you can use this code to make a gif
% convert -loop 0 -delay 1 1*.png movie.gif

%% Save the cams to a text file


%%%%%%%%%%% Homework assignment
% write the data into files(left and right cam)
% make up z data point
z0 = zeros(length(x_cam_left),1);
dlmwrite('x_y_data_leftCAM.txt',[x_cam_left y_cam_left z0],',');
dlmwrite('x_y_data_rightCAM.txt',[x_cam_right y_cam_right z0],',');

