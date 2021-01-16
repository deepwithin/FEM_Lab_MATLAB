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

%% meshing of the geometry
n_members = 3; % no. of truss members
elements_per_member = 2;
ne = n_members * elements_per_member; % total no of elements
%nne = 2; % no. of noded element
le = L/elements_per_member; % length of each element (if all equal)
nn = ne; % total no of nodes = no. of elemets in truss
dofn = 2; % no. of degrees of freedom per node
dofe = 4; % dof per element (2 noded element)
tdof = dofn * nn; % total degrees of freedom
ith_member = 1; %member check for diff- E, A, le


%% connectivity matrix generator (relation of local no with global numbers)
conn = zeros(ne,dofe);
fill=1;
for i4 = 1:ne
    for i5 = 1:dofe
        conn(i4,i5) = fill;
        fill = fill + 1;
        if fill > tdof && i4 == ne
            fill = 1;
        end
    end
    fill = fill - dofn;
end

%% initializing vectors
KG = zeros(tdof,tdof); % initializing global stiffness matrix
FGU = zeros(tdof,1); % initalization of uniformly distributed load vector
FGC = zeros(tdof,1); % initialization of concentrated load vector

%% populating global stiffness matrix
for i = 1:ne
    if i > elements_per_member * ith_member
        ith_member = ith_member + 1;
    end
    c = cosd(element_angle(ith_member));
    s = sind(element_angle(ith_member));
    Ke = E(ith_member)*A(ith_member)/le(ith_member) * [c*c c*s -c*c -c*s;...
                                                       c*s s*s -c*s -s*s;...
                                                       -c*c -c*s c*c c*s;...
                                                       -c*s -s*s c*s s*s]
    Fe = zeros(dofe,1);
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j), conn(i,k)) =  KG(conn(i,j), conn(i,k)) + Ke(j,k); 
        end
        FGU(conn(i,j),1) =  FGU(conn(i,j),1)+ Fe(j,1);
    end
end

%% load and reaction vector
FGC(end-1) = Px; % applying load values
FGC(end) = Py;
FG = FGU + FGC; % global force vector
KG_nbc = KG; % K global with no boundary conditions

%% applying boundary conditions
for i2 = 1:((elements_per_member+1) * dofn)
   KG(i2,:) = 0;
   KG(:,i2) = 0;
   KG(i2,i2) = 1;
   
end
%KG(7,7)=1;
%% solving the matrix equation
UG = linsolve(KG,FG)

reaction_vector = KG_nbc * UG - FG

%% Convergence study


%% Graph plotting
  



%% Function of entire process for convergence
