% clc;
% clear all;
% close all;
x = 1;
y = 1;
dx = 0.025;
dy = 0.025;

Nx = x/dx + 1;
Ny = y/dy + 1;

j1 = 20;
j2 = 25;

j3 = 5;
j4 = 10;

psi_1 = -0;
%u0 = 1;

psi = zeros(Ny,Nx);
w   = zeros(Ny,Nx);
u   = zeros(Ny,Nx)+u0;
v   = zeros(Ny,Nx);

Re = 500;
gamma = u0/Re;
%alpha = 1.8;
err = 10;
err1 = err;
iter = 1;


dt = 0.02;