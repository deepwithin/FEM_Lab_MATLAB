%% Static analysis of a Bar
clc
clear all
close all

%Geometrical Parameters
E = 200*10^9;
L = 50e-2;
A = 25e-4;
P = 2000;

%Meshing
ne = 5;
nne = 2;
le = L/ne;
nn = ne + 1;
dofn = 1;
dofe = 2;
tdof = dofn*nn;
conn = [1 2;2 3;3 4;4 5;5 6]; %connectivity Matrix

KG = zeros(tdof,tdof);
%FG = zeros(tdof,1);
FGU = zeros(tdof,1);
FGC = zeros(tdof,1);
for i= 1:ne
    Ke = E*A/le*[1 -1;-1 1];
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