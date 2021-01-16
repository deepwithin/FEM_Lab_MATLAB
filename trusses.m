%%%%%%%%%%%%%%%%%%%%%%%%
%%_TRUSS FEM ANALYSIS_%%
%%%%%%%%%%%%%%%%%%%%%%%%

clc %clears the screen
clear all %closes all the figure windows
close all %closes all the fig windows

%% geometrical parameters and material properties
E = 10^9*[200 200 200]; % young's modulus in Pa
L = 10^-2*[100 100 100*sqrt(2)]; % entire length of beam
A = 10^-4*[30 30 30]; % area of cross sections
Px = 100*10^3; % conc. load value in x dirn
Py = -200*10^3; % conc. load value in y dirn

element_angle = [0 90 225]; % element orientation angles

X = [0 L(1) L(1) 0];
Y = [0 0 L(1) 0];

%% meshing of the geometry
ne = 3; % total no of elements
nne = 2; % no. of noded element
nn = ne; % total no of nodes = no. of elemets in truss
dofn = 2; % no. of degrees of freedom per node
dofe = 4; % dof per element (2 noded element)
tdof = dofn * nn; % total degrees of freedom
conn  = [3 4 5 6;...
         5 6 1 2;...
         3 4 1 2]; % connectivity matrix (relation of local no with global numbers)

%% initializing vectors
KG = zeros(tdof,tdof); % initializing global stiffness matrix
FGU = zeros(tdof,1); % initalization of uniformly distributed load vector
FGC = zeros(tdof,1); % initialization of concentrated load vector

%% populating global stiffness matrix
for i = 1:ne
    c = cosd(element_angle(i));
    s = sind(element_angle(i));
    Ke = E(i)*A(i)/L(i) * [c*c c*s -c*c -c*s;...
                           c*s s*s -c*s -s*s;...
                           -c*c -c*s c*c c*s;...
                           -c*s -s*s c*s s*s];
    Fe = zeros(dofe,1);
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j), conn(i,k)) =  KG(conn(i,j), conn(i,k)) + Ke(j,k); 
        end
        FGU(conn(i,j),1) =  FGU(conn(i,j),1)+ Fe(j,1);
    end
end

%% load and reaction vector
FGC(5) = Px; % applying load values
FGC(6) = Py;

FG = FGU + FGC; % global force vector
FG_nbc = FG; % F global with no boundary conditions
KG_nbc = KG; % K global with no boundary conditions

%% applying boundary conditions
for i2 = 1:4
   KG(i2,:) = 0;
   KG(:,i2) = 0;
   KG(i2,i2) = 1;
   
end

%% solving the matrix equation
UG = linsolve(KG,FG)

R = KG_nbc * UG - FG_nbc

%% Convergence study


%% Visualization
hold on

X_new = X;
Y_new = Y;
plot(0,0,'r*');
plot(X,Y,'g-o');

amplification = 100;
i4=1;
for i3 = [1,3,5]
    X_new(i4) = X(i4) + amplification*UG(i3);
    i4 = i4+1;
end
i4=1;
for i3 = [2,4,6]
    Y_new(i4) = Y(i4) + amplification*UG(i3);
    i4 = i4+1;
end
plot(X_new,Y_new,'b--o');
hold off
