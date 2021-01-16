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
E = 200*10^9; % young's modulus in Pa
v = 0.3; % poisson's ratio
h = 2 * 10^-2; % plate thickness (element wise)


global X
X = [0 1 2  0 1 2  0 1 2];
global Y
Y = [0 0 0  1 1 1  2 2 2];

% for dynamic generation of conn, L, angles, etc
global econn
econn = [1 2 4 5;...
         2 3 5 6;...
         4 5 7 8;...
         5 6 8 9]; %element connectivity (anti-clock wise)


%% meshing of the geometry
ne = 4; % total no of elements
nne = 4; % no. of noded element

nn = 9; % total no of nodes = no. of elemets in truss
dofn = 2; % no. of degrees of freedom per node
dofe = 8; % dof per element (2 noded element)
tdof = dofn * nn; % total degrees of freedom
conn  = [1 2  3 4  7 8  9 10;...
         3 4  5 6  9 10  11 12;...
         7 8  9 10  13 14  15 16;...
         9 10  11 12  15 16  17 18]; % connectivity matrix (relation of local no with global numbers)

%% initializing vectors
KG = zeros(tdof,tdof); % initializing global stiffness matrix
FGU = zeros(tdof,1); % initalization of uniformly distributed load vector
FGC = zeros(tdof,1); % initialization of concentrated load vector

%% populating global stiffness matrix
for i = 1:ne
    
    he = h;
    Ke = zeros(dofe,dofe);
    Fe = zeros(dofe,1);
    
    [Xi_vec, Eta_vec, Weight] = gauss_legendre_calc(2,2);
    
    x1 = X(econn(i,1));
    x2 = X(econn(i,2));
    x3 = X(econn(i,3));
    x4 = X(econn(i,4));
  
    y1 = Y(econn(i,1));
    y2 = Y(econn(i,2));
    y3 = Y(econn(i,3));
    y4 = Y(econn(i,4));
    
    D = E/(1-v^2) * [1  v  0;...
                     v  1  0;...
                     0  0  (1-v)/2]; %constitutive relations matrix
    
    for i6 = 1:4

        Psi1_Xi_value = (1/4)*[(Eta_vec(i6)-1),0,(1-Eta_vec(i6)),0,-(1+Eta_vec(i6)),0,(1+Eta_vec(i6)),0];
        Psi1_Eta_value = (1/4)*[-(1-Xi_vec(i6)),0,-(1+Xi_vec(i6)),0, (1-Xi_vec(i6)),0,(1+Xi_vec(i6)),0];
        Psi2_Xi_value = (1/4)*[0,(Eta_vec(i6)-1),0,(1-Eta_vec(i6)),0,-(1+Eta_vec(i6)),0,(1+Eta_vec(i6))];
        Psi2_Eta_value = (1/4)*[0,-(1-Xi_vec(i6)),0,-(1+Xi_vec(i6)),0, (1-Xi_vec(i6)),0,(1+Xi_vec(i6))];


        Pe = [x1; y1; x2; y2; x3; y3; x4; y4];
        
        J = [];
        
        J(1,1) = Psi1_Xi_value*Pe;
        J(1,2) = Psi2_Xi_value*Pe;
        J(2,1) = Psi1_Eta_value*Pe;
        J(2,2) = Psi2_Eta_value*Pe;

        J_det = det(J);

        A = (1/J_det)*[J(2,2),     -J(1,2),       0,             0; ...
                       0,          0,             -J(2,1),       J(1,1);...
                       -J(2,1),    J(1,1),        J(2,2),        -J(1,2)];

        G = [Psi1_Xi_value;...
             Psi1_Eta_value;...
             Psi2_Xi_value;...
             Psi2_Eta_value];

        B = A*G;

        Ke_nodal = he*( (transpose(B))*D*B )* J_det* Weight(i6);
        
        Ke = Ke + Ke_nodal %global elemental stiffness matrix

    end
                 
    
    le = 0;
%     Fe = h*le/6 * [2*Tx1 + Tx2;...
%                       2*Ty1 + Ty2;...
%                       Tx1 + 2*Tx2;...
%                       Ty1 + 2*Ty2;...
%                       0;...
%                       0];%elemental load vector
    
    for j = 1:dofe
        for k = 1:dofe
            KG(conn(i,j), conn(i,k)) =  KG(conn(i,j), conn(i,k)) + Ke(j,k); 
        end
        FGU(conn(i,j),1) =  FGU(conn(i,j),1)+ Fe(j,1);
    end
end

%% load and reaction vector
FGC(18) = -30000;

FG = FGU + FGC; % global force vector
FG_nbc = FG; % F global with no boundary conditions
KG_nbc = KG; % K global with no boundary conditions

%% applying boundary conditions
for i2 = [1,2,7,8,13,14]
   KG(i2,:) = 0;
   KG(:,i2) = 0;
   KG(i2,i2) = 1;
   
   FG(i2,1) = 0;
end

%% solving the matrix equation
UG = linsolve(KG,FG)
fprintf('\n');
R = KG_nbc * UG - FG_nbc

% disp('displacement at A, x value');
% UG(5)
% disp('displacement at A, y value');
% UG(6)
% 
% disp('displacement at B, x value');
% UG(3)
% disp('displacement at B, y value');
% UG(4)
% 
% 
% disp('reaction at C, x value');
% R(1)
% disp('reaction at C, y value');
% R(2)
% 
% disp('reaction at D, x value');
% R(7)
% disp('reaction at D, y value');
% R(8)

% KG

%% Convergence study


%% Visualization
hold on

X_new = X;
Y_new = Y;
plot(0,0,'r*');
plot(X,Y,'g-o');

amplification = 500;
i4=1;
for i3 = 1:2:17
    X_new(i4) = X(i4) + amplification*UG(i3);
    i4 = i4+1;
end
i4=1;
for i3 = 2:2:18
    Y_new(i4) = Y(i4) + amplification*UG(i3);
    i4 = i4+1;
end
plot(X_new,Y_new,'b--o');
hold off

%% Functions
function [xi,eta,weight] = gauss_legendre_calc(x_point,y_point)

ngx = x_point;
ngy = y_point;

GP = [0,                 0,                 0,                 0;...
      -1/sqrt(3),        1/sqrt(3),         0,                 0;...
      -sqrt(3/5),        0,                 sqrt(3/5),         0;...
      -0.8611363116,    -0.3399810435,      0.3399810435,      0.8611363116]; %gauss points
  
WP = [2,                   0,                  0,                0;...
      1,                   1,                  0,                0;...
      5/9,                 8/9,                5/9,              0;...
      0.3478548451,        0.6521451548,       0.3478548451,     0.6521451548]; %weight points
  
Gauss_point = GP(ngx,:);

W_x = WP(ngx,:);
W_y = WP(ngy,:);


Xi_r = [];
Eta_r = [];
Weight_r = [];

for ix = 1:ngx
    for iy = 1:ngy
        Weight_r(ix,iy) = W_x(ix) * W_y(iy);
   
        Xi_r(ix,iy)= Gauss_point(ix);
        Eta_r(ix,iy) = Gauss_point(iy);
   
    end
end

Xi_r_reshape = reshape(Xi_r,1,ngx*ngx);
Eta_r_reshape = reshape(Eta_r,1,ngy*ngy);
Weight_r_reshape = reshape(Weight_r,1,ngx*ngy);

xi = Xi_r_reshape;
eta = Eta_r_reshape;
weight = Weight_r_reshape;

end