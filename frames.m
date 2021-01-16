%%%%%%%%%%%%%%%%%%%%%%%%%
%%_FRAMES FEM ANALYSIS_%%
%%%%%%%%%%%%%%%%%%%%%%%%%

clc %clears the screen
clear all %closes all the figure windows
close all %closes all the fig windows
format
%format shortEng
%format shortG
format compact
%format loose

%% geometrical parameters and material properties
E = 10^9 * 200; % young's modulus in Pa
L = 10^-2 * [40 25 40]; % entire length of beam
A = 4 * 10^4 * 10^-6; % area of cross sections
I = 20 * 10^6 * 10^-12; % moment of area

q0 = [0 -10000 0]; %uniform loading

element_angle = [90 0 270];

X = [0 0 L(2) L(2)];
Y = [0 L(1) L(3) 0];

% for dynamic generation of conn, L, angles, etc
econn = [1 2;...
         2 3;...
         3 4]; %element connectivity

angle = zeros(1,3);
for i1 = 1:3
    angle(i1) = atand(( Y(econn(i1,2)) - Y(econn(i1,1)) )/( X(econn(i1,2)) - X(econn(i1,1)) ));
end

%% meshing of the geometry
ne = 3; % total no of elements
nne = 2; % no. of noded element
Le = L; % length of each element
nn = ne+1; % total no of nodes = no. of elemets in truss
dofn = 3; % no. of degrees of freedom per node
dofe = 6; % dof per element (2 noded element)
tdof = dofn * nn; % total degrees of freedom
conn  = [1 2 3  4 5 6;...
         4 5 6  7 8 9;...
         7 8 9  10 11 12]; % connectivity matrix (relation of local no with global numbers)

%% initializing vectors
KG = zeros(tdof,tdof); % initializing global stiffness matrix
FGU = zeros(tdof,1); % initalization of uniformly distributed load vector
FGC = zeros(tdof,1); % initialization of concentrated load vector

%% populating global stiffness matrix
for i = 1:ne
    c = cosd(element_angle(i));%cos of element angle in degrees
    s = sind(element_angle(i));%sin of element angle in degrees
    le = Le(i);
    ax = E*A/le; %axial component
    tr = E*I/(le^3); %transverse component
    
    Ke = [ax       0           0             -ax       0            0;...
          0        12*tr       6*le*tr       0         -12*tr       6*le*tr;...
          0        6*le*tr     4*le^2*tr     0         -6*le*tr     2*le^2*tr;...
          -ax      0           0             ax        0            0;...
          0        -12*tr      -6*le*tr      0         12*tr        -6*le*tr;...
          0        6*le*tr     2*le^2*tr     0         -6*le*tr     4*le^2*tr];%elemental local stiffness matrix
    
    T = [c  s  0  0  0  0;...
         -s c  0  0  0  0;...
         0  0  1  0  0  0;...
         0  0  0  c  s  0;...
         0  0  0 -s  c  0;...
         0  0  0  0  0  1]; %global to local transformation matrix
    
    Ke = T.' * Ke * T; %global elemental stiffness matrix
    
     %Fe = zeros(dofe,1);
    Fe = q0(i)*le/2 * [0; 1; le/6; 0; 1; -le/6];%elemental load vector
    Fe = T.' * Fe; %global elemental load vector
    
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j), conn(i,k)) =  KG(conn(i,j), conn(i,k)) + Ke(j,k); 
        end
        FGU(conn(i,j),1) =  FGU(conn(i,j),1)+ Fe(j,1);
    end
end

%% load and reaction vector
%FGC(4) = 0;

FG = FGU + FGC; % global force vector
FG_nbc = FG; % F global with no boundary conditions
KG_nbc = KG; % K global with no boundary conditions

%% applying boundary conditions
for i2 = [1,2,3,10,11,12]
   KG(i2,:) = 0;
   KG(:,i2) = 0;
   KG(i2,i2) = 1;
   
   FG(i2,1) = 0;
end

%% solving the matrix equation
UG = linsolve(KG,FG)
fprintf('\n');
R = KG_nbc * UG - FG_nbc

%% Convergence study


%% Visualization
hold on

X_new = X;
Y_new = Y;
plot(0,0,'r*');
plot(X,Y,'g-o');

amplification = 100000;
i4=1;
for i3 = [1,4,7,10]
    X_new(i4) = X(i4) + amplification*UG(i3);
    i4 = i4+1;
end
i4=1;
for i3 = [2,5,8,11]
    Y_new(i4) = Y(i4) + amplification*UG(i3);
    i4 = i4+1;
end
plot(X_new,Y_new,'b--o');
hold off
