close all
clear all
clc

Ts=0.001;       %sample time
t=0:Ts:20;   %time vector

%starting point     qi=[1+alpha 1 pi/4] -> P38000219 -> alpha=9
x_i = 10;       
y_i = 1;
theta_i = pi/4 ;

%final point
x_f= 0;
y_f = 0;
theta_f = 0; 

%k values (chosen the same as professor Siciliano's book)
k1=1;
k2=2.5;
k3=3;


simOut =sim('es5_RL',t);    %start the simulation directly from the script
x_k=simOut.x_k;
y_k=simOut.y_k;
theta_k=simOut.theta_k;
x=simOut.x;
y=simOut.y;
theta=simOut.theta;

v=simOut.vel;
w=simOut.omega;


%Plots from Runge-Kutta model approximation
figure
plot(t,theta_k);
xlabel('time');
ylabel('theta');
title('Robot orientation');

figure
plot(x_k,y_k)
xlabel('x');
ylabel('y');
title('R-K Approximated Path Trajectory');

figure
plot(v);
xlabel('time');
ylabel('v(t)');
title('Linear velocity');

figure
plot(w);
xlabel('time');
ylabel('w(t)');
title('Angular velocity');

%Plot from the real unicycle
figure
plot(x,y);
xlabel('x');
ylabel('y');
title('Robot Path Trajectory');
