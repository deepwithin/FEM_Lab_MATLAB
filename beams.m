%%%%%%%%%%%%%%%%%%%%%%%%
%%_BEAMS FEM ANALYSIS_%%
%%%%%%%%%%%%%%%%%%%%%%%%

clc %clears the screen
clear all %closes all the figure windows
close all %closes all the fig windows
format
%format shortEng
%format shortG
format compact
%format loose

%% geometrical parameters and material properties
E = 10^9*200; % young's modulus in Pa
L = 2; % entire length of beam
I = 4*10^6*10^-12; % moment of area
Px = 0; % conc. load value in x dirn
Py = 0; % conc. load value in y dirn
q0 = [0 -12000]; %uniform loading

X = [0 1 2];
slope = [0 0 0];
Y = [0 0 0];

%% meshing of the geometry
ne = 2; % total no of elements
nne = 2; % no. of noded element
le = L/ne; % length of each element
nn = ne+1; % total no of nodes = no. of elemets in truss
dofn = 2; % no. of degrees of freedom per node
dofe = 4; % dof per element (2 noded element)
tdof = dofn * nn; % total degrees of freedom
conn  = [1 2 3 4;...
         3 4 5 6;...
         1 2 5 6]; % connectivity matrix (relation of local no with global numbers)

%% initializing vectors
KG = zeros(tdof,tdof); % initializing global stiffness matrix
FGU = zeros(tdof,1); % initalization of uniformly distributed load vector
FGC = zeros(tdof,1); % initialization of concentrated load vector

%% populating global stiffness matrix
for i = 1:ne
    
    Ke = E*I/(le^3) * [12      6*le      -12      6*le;...
                       6*le    4*le^2    -6*le    2*le^2;...
                       -12     -6*le^2   12       -6*le;...
                       6*le    2*le^2    -6*le    4*le^2]; %elemental stiffness matrix
    %Fe = zeros(dofe,1);
    Fe = q0(i)*le/2 * [1; 1/6; 1; -1/6];
    
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j), conn(i,k)) =  KG(conn(i,j), conn(i,k)) + Ke(j,k); 
        end
        FGU(conn(i,j),1) =  FGU(conn(i,j),1)+ Fe(j,1);
    end
end

%% load and reaction vector

FG = FGU + FGC; % global force vector
FG_nbc = FG; % F global with no boundary conditions
KG_nbc = KG; % K global with no boundary conditions

%% applying boundary conditions
for i2 = [1,2,3,5]
   KG(i2,:) = 0;
   KG(:,i2) = 0;
   KG(i2,i2) = 1;
   
   FG(i2,1) = 0;
end

%% solving the matrix equation
UG = linsolve(KG,FG)

R = KG_nbc * UG - FG_nbc

%% Convergence study


%% Visualization
hold on

X_new = X;
slope_new = slope;
Y_new = Y;
plot(0,0,'r*');
plot(X,Y,'g-o');

amplification = 1;
i4=1;
for i3 = [1,3,5]
    Y_new(i4) = Y(i4) + amplification*UG(i3);
    i4 = i4+1;
end
i4=1;
for i3 = [2,4,6]
    slope_new(i4) = slope(i4) + amplification*UG(i3);
    i4 = i4+1;
end
% for i3 = 1 : length(X)
%     Y_new(i3) = X(i3)*slope_new(i3);
% end
plot(X_new,Y_new,'b--o');
hold off
