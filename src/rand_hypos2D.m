function [rxp,ryp] = rand_hypos2D(L,zerr_av,nhypos,strike,rt)

%  construct set of random hypocenters on the line with an additional
%  random error.  The line is specified by -L/2 < x < L/2 and
%  strike direction.  nhypos is the number of hypocenters desired.

%  zerr_av is the standard deviation of z error
%  strike, dip, and rake are in the geographical coordinate system
%  (degrees)

%  rt is a constant translation vector of the origin.
% clear all; close all; clc;
% 
% L= 100; zerr_av=10;rt = [0,0];
% nhypos = 200; strike=360;


con=pi/180.;
strike=strike.*con;

%close all;

n=nhypos/4;
ry=[ (randperm(n) + randn(1,n))  (randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n))] ;

% scale to L
ry=ry.*L./(2.0.*n);

% now construct an error in z
rx0=randn(1,nhypos);
rx = (rx0/max(rx0))*zerr_av/2;

% rotate into strike direction
R=[rx ; ry];
Dstrike=[ cos(strike)  sin(strike); -sin(strike) cos(strike)];

Rstrike=Dstrike*R;

rxp(1:nhypos)=Rstrike(1,1:nhypos);
ryp(1:nhypos)=Rstrike(2,1:nhypos);

%  translate the final rotated collection of points
rxp=rxp+rt(1);
ryp=ryp+rt(2);
