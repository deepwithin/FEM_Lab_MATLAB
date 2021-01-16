%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%_PLANE STRESS FEM ANALYSIS_%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc %clears the screen
clear all %closes all the figure windows
close all %closes all the fig windows
format
%format shortEng
%format shortG
format compact
%format loose

%% geometrical parameters and material properties
E = 200*10^9 * [1 1 1 1 1 1 1 1 1]; % young's modulus in Pa
v = 0.25; % poisson's ratio
h = [0.01 0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01]; % plate thickness (element wise)

theta = atand(60/30);

Tx1 = 10^6 * [0 0 0 1*(-cosd(theta)) 2*(-cosd(theta)) 0 0 0];
Tx2 = 10^6 * [0 0 0 2*(-cosd(theta)) 3*(-cosd(theta)) 0 0 0];
Ty1 = 10^6 * [0 0 0 1*(-sind(theta)) 2*(-sind(theta)) 0 0 0];
Ty2 = 10^6 * [0 0 0 2*(-sind(theta)) 3*(-sind(theta)) 0 0 0];

global X
X = [0 70 100 0 50 85 0 25 70];
global Y
Y = [20 20 20 40 40 40 60 60 60];

% for dynamic generation of conn, L, angles, etc
global econn
econn = [1 5 4;...
         2 5 1;...
         3 6 2;...
         6 9 5;...
         5 9 8;...
         5 8 4;...
         4 8 7]; %element connectivity (anti-clock wise)


%% meshing of the geometry
ne = 2; % total no of elements
nne = 3; % no. of noded element

nn = 4; % total no of nodes = no. of elemets in truss
dofn = 2; % no. of degrees of freedom per node
dofe = 6; % dof per element (2 noded element)
tdof = dofn * nn; % total degrees of freedom
conn  = [1 2  9 10  7 8;...
         3 4  9 10  1 2;...
         3 4  11 12  9 10;...
         ]; % connectivity matrix (relation of local no with global numbers)

%% initializing vectors
KG = zeros(tdof,tdof); % initializing global stiffness matrix
FGU = zeros(tdof,1); % initalization of uniformly distributed load vector
FGC = zeros(tdof,1); % initialization of concentrated load vector

%% populating global stiffness matrix
for i = 1:ne
    
    s1 = sqrt( x_ab(1,2,i)^2 + y_ab(1,2,i)^2 );
    s2 = sqrt( x_ab(2,3,i)^2 + y_ab(2,3,i)^2 );
    s3 = sqrt( x_ab(3,1,i)^2 + y_ab(3,1,i)^2 );
    s = (s1+s2+s3)/2;
    A = sqrt(s*(s-s1)*(s-s2)*(s-s3)); %area by heron's formula
    
    J = [x_ab(1,3,i)      y_ab(1,3,i);...
         x_ab(2,3,i)      y_ab(2,3,i)]; %elemental Jacobian matrix
    
    %strain displacement relation matrix
    B = 1/det(J)*[y_ab(2,3,i)  0            y_ab(3,1,i)   0            y_ab(1,2,i)  0;...
                  0            x_ab(3,2,i)  0             x_ab(1,3,i)  0            x_ab(2,1,i);...
                  x_ab(3,2,i)  y_ab(2,3,i)  x_ab(1,3,i)   y_ab(3,1,i)  x_ab(2,1,i)  y_ab(1,2,i)];
    
    D = E(i)/(1-v^2) * [1  v  0;...
                        v  1  0;...
                        0  0  (1-v)/2]; %constitutive relations matrix
    
    Ke = A*h(i) * B.' * D * B; %global elemental stiffness matrix
    
%     Fe = zeros(dofe,1);
    le = 0;
    Fe = h(i)*le/6 * [2*Tx1 + Tx2;...
                      2*Ty1 + Ty2;...
                      Tx1 + 2*Tx2;...
                      Ty1 + 2*Ty2;...
                      0;...
                      0];%elemental load vector
    
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j), conn(i,k)) =  KG(conn(i,j), conn(i,k)) + Ke(j,k); 
        end
        FGU(conn(i,j),1) =  FGU(conn(i,j),1)+ Fe(j,1);
    end
end

%% load and reaction vector
FGC(4) = -1000;

FG = FGU + FGC; % global force vector
FG_nbc = FG; % F global with no boundary conditions
KG_nbc = KG; % K global with no boundary conditions

%% applying boundary conditions
for i2 = [2,5,6,7,8]
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

amplification = 1000;
i4=1;
for i3 = [1,3,5,7]
    X_new(i4) = X(i4) + amplification*UG(i3);
    i4 = i4+1;
end
i4=1;
for i3 = [2,4,6,8]
    Y_new(i4) = Y(i4) + amplification*UG(i3);
    i4 = i4+1;
end
plot(X_new,Y_new,'b--o');
hold off

%% Functions
function outx = x_ab(a,b,i)
    global X
    global econn
    outx = X(econn(i,a)) - X(econn(i,b));
end

function outy = y_ab(a,b,i)
    global Y
    global econn
    outy = Y(econn(i,a)) - Y(econn(i,b));
end