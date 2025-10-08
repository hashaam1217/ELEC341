% Date: 08/10/2025
% Name: Hashaam Zafar
% Student Number: 10078020

% Setup
clear all; clc;
SN = 10078020;

A = 11;
B = 10; 
C = 10; 
D = 17; 
E = 18;
F = 10; 
G = 12; 
H = 10; 

s = tf('s');

% Figure 1 and 2
G1 = 1 / (s + A);
G2 = 5 / (s + B);
G3 = 5 / (s + C);
G4 = 5 / (s + D);
G5 = 1 / (s + E);
H1 = 100 / (s + F); 
H2 = 500 / (s + G); 

% Figure 3
Jf = A / 25; 
Bf = B / 20; 
Jg = C * 50; 
Bg = D / 30; 
n = E; 
