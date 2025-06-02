clear
clc
close all

ha = 6000;
hp = 1200;
u = 398600;
R = 6378;

ra = ha + R;
rp = hp + R;

e = (ra - rp)/(ra + rp);
a = (ra + rp)/2;
P = a * (1 - e^2);
P = 14850;
e = 0.1;

r = @(theta) P / (1 + e*cos(theta));
vr = @(theta) sqrt(u/P) * e * sin(theta);
v0 = @(theta) sqrt(u/P) * (1 + e*cos(theta));

theta = 180 * pi/180;
omega = 45 * pi/180;
w = 30 * pi/180;
i = 15 * pi/180;

r1 = r(theta);
vr1 = vr(theta);
v01 = v0(theta);

r1 = [r1; 0; 0];
v1 = [vr1; v01; 0];

R_perifocale = @(theta) [cos(theta), sin(theta), 0;
                         -sin(theta), cos(theta), 0;
                         0, 0 , 1];

R1 = @(omega) [cos(omega), sin(omega), 0;
               -sin(omega), cos(omega), 0;
               0, 0, 1];
R2 = @(i) [1, 0, 0;
           0, cos(i), sin(i);
           0, -sin(i), cos(i)];
R3 = @(w) [cos(w), sin(w), 0;
           -sin(w), cos(w), 0;
           0, 0, 1];
R = @(omega, i, w) R1(omega)*R2(i)*R3(w);


r1_perifocale = R_perifocale(theta) * r1;
v1_perifocale = R_perifocale(theta) * v1;

r1_geocentrico = R(omega, i, w) * r1_perifocale;
v1_geocentrico = R(omega, i, w) * v1_perifocale;

