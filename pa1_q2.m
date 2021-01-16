%% Static analysis of a Bar
clc
clear all
close all

%Geometrical Parameters
E = 200*10^9;
L = 50e-2;
dia_P = 0.15;
dia_Q = 0.07;
P = 5000;

%Meshing
ne = 5;
nne = 2;
le = L/ne;
nn = ne + 1;
dofn = 1;
dofe = 2;
tdof = dofn*nn;
x=0;

conn = [1 2;2 3;3 4;4 5;5 6];

KG = zeros(tdof,tdof);
%FG = zeros(tdof,1);
FGU = zeros(tdof,1);
FGC = zeros(tdof,1);
for i= 1:ne
    Ke = E*(pi/4*(dia_of_element(i, le, dia_P, dia_Q, L))^2)/le*[1 -1;-1 1];
    x=x+1;
    Fe = zeros(2,1);
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j),conn(i,k)) = KG(conn(i,j),conn(i,k)) + Ke(j,k);
        end
        FGU(conn(i,j),1) = FGU(conn(i,j),1) + Fe(j,1);
    end
end

FGC(end) = 2000; %concentrated load vector

FG = FGU + FGC;


% Application of boundary conditions
KG(1,:) = 0;
KG(:,1) = 0;
KG (1,1) = 1;
FG(1,1) = 0;

UG = linsolve(KG,FG)


% function definitions

function dia = dia_at_x(x, dia_P, dia_Q, L)
    width = 15-(dia_P-dia_Q)/L*x;
end

function dia = dia_of_element(element_number, le, dia_P, dia_Q, L)
    dia = (dia_at_x((element_number-1)*le, dia_P, dia_Q, L)+dia_at_x(element_number*le, dia_P, dia_Q, L))/2;
end