close all; clear all; clc;
%%% Solver for the inhomogeneous model problem in Mani and Park (2021) %%%
% staggered mesh using 2nd order finite difference
% scalar c is stored at cells
% velocities u1 and u2 are stored at faces


%% Parameters
DM1 = 0.05;
DM2 = 1;

% number of mesh points
N1 = 100; %40;   
N2 = N1/2; %10;    

Nvx = 1;
Nvy = 1;

if mod(Nvx,2) == 0
    phaseskipX = pi
else
    phaseskipX = 0
end

%% Mesh Generation (staggered mesh)
dx1 = 2*pi/N1;
dx2 = 2*pi/N2;
% cells
x1 = -pi+dx1/2:dx1:pi-dx1/2;  
x2 = dx2/2:dx2:2*pi-dx2/2;
% faces
x1_f = -pi:dx1:pi; 
x2_f = 0:dx2:2*pi; 

%% Velocity field
% create u1
display('get u1')
u1 = zeros(N1+1,N2);
for i=1:N1+1
    for j=1:N2
%         u1(i,j) = (1+cos(x1_f(i)))*cos(x2(j));
        u1(i,j) = (1+cos(phaseskipX+Nvx*x1_f(i)))*cos(Nvy*x2(j));
    end
end

% create u2
display('get u2')
u2 = zeros(N1,N2+1);
% u2 that enforces incompressibility
for i=1:N1
    for j=2:N2
        u2(i,j) = u2(i,j-1) - dx2/dx1*(u1(i+1,j-1)-u1(i,j-1)); 
    end
end

disp('reshape u1 u2')
%% reshape u1 and u2 
u1 = u1(:);
u1 = spdiags(u1,0,size(u1,1),size(u1,1));

u2 = u2(:);
u2 = spdiags(u2,0,size(u2,1),size(u2,1));


%% Operators
disp('operators')
% interpolation matrix from center to edge in x1
% size(diag(0.5*ones(N1,1),-1))
% size(diag(0.5*ones(N1+1,1),0))
% 
% size(spdiags(0.5*ones(N1,1),-1,N1,N1))
% size(spdiags(0.5*ones(N1+1,1),0,N1,N1))

interp1 = spdiags(0.5*ones(N1,1),-1,N1+1,N1+1) + spdiags(0.5*ones(N1+1,1),0,N1+1,N1+1);
interp1 = interp1(:,1:N1);

% add dirichlet bcs  
interp1(1,1) = 0;
interp1(N1+1,N1) = 0;
interpx1 = kron(speye(N2), interp1);


% create d/dx1 from edge to center
ddx1_ec = spdiags(ones(N1+1,1),1,N1+1,N1+1) + spdiags(-ones(N1+1,1),0,N1+1,N1+1);
ddx1_ec = (1/dx1)*ddx1_ec(1:N1,:);
ddx1 = kron(speye(N2), ddx1_ec);


% interpolation matrix from center to edge in x2
interpx2 = kron(spdiags(0.5*ones(N2+1,1),-1,N2+1,N2+1),speye(N1)) + ...
    kron(spdiags(0.5*ones(N2+1,1),0,N2+1,N2+1),speye(N1));
interpx2 = interpx2(:,1:N1*N2);
% add dc/dx2=0 bcs
interpx2(1:N1,1:N1) = speye(N1);
interpx2(N1*N2+1:N1*(N2+1),N1*(N2-1)+1:N1*N2) = speye(N1);

% create d/dx2 from edge to center
ddx2 = kron(spdiags(ones(N2+1,1),0,N2+1,N2+1),-speye(N1)) + ...
    kron(spdiags(ones(N2+1,1),1,N2+1,N2+1),speye(N1));
ddx2 = (1/dx2)*ddx2(1:N1*N2,:);

% create d^2/(dx1)^2
d2dx12_bar = spdiags(ones(N1,1),-1,N1,N1) + spdiags(-2*ones(N1,1),0,N1,N1) + spdiags(ones(N1,1),1,N1,N1);
% add dirichlet BCs, 2nd order interpolation
d2dx12_bar(1,1) = d2dx12_bar(1,1)-2;
d2dx12_bar(1,2) = d2dx12_bar(1,2)+1/3;
d2dx12_bar(N1,N1) = d2dx12_bar(N1,N1)-2; 
d2dx12_bar(N1,N1-1) = d2dx12_bar(N1,N1-1)+1/3; 
d2dx12_bar = 1/dx1^2*d2dx12_bar;
d2dx12 = kron(speye(N2), d2dx12_bar); 

% create d^2/(dx2)^2
% diagonal elements
Dia = kron(speye(N2), spdiags(-2*ones(N1,1),0,N1,N1)); 
% offdiagonal elements
ODia = kron(spdiags(ones(N2,1),-1,N2,N2)+spdiags(ones(N2,1),1,N2,N2), speye(N1));
% add dc/dx2=0 BCs
Dia(1:N1,1:N1) = -speye(N1); 
Dia(N1*(N2-1)+1:N1*N2,N1*(N2-1)+1:N1*N2) = -speye(N1); 
d2dx22 = 1/dx2^2*(Dia+ODia);

%% L matrix
display('assemble L')
L = ddx1*(u1*interpx1) + ddx2*(u2*interpx2) - DM1*d2dx12 - DM2*d2dx22; 

%% P Matrix   
% averaged in x2 direction
display('P matrix')
n = N2;
Pc = kron(ones(1,N2), 1/n*speye(N1));

%% E Matrix
display('E matrix')
Ec = n*Pc';

% error()
%% L_bar Matrix
% L_inv = inv(L);
disp('getting lec');
LEC = L\Ec;
disp('getting lbar');
L_bar = inv(Pc*LEC);
disp('got lbar');

%% Plot
figure()
pcolor(log10(abs(L_bar)))
axis ij
caxis([-2 2])
colorbar
shading flat
title('$\bar{\mathcal{L}}$','Interpreter','Latex')

index = size(L_bar,1)/2;
figure()
plot(log10(abs(L_bar(index,:))))
% error()


%% Eddy-Diffusivity 
% ddx1 cell center to edge
ddx1_ce = spdiags(-ones(N1+1,1),-1,N1+1,N1+1)+spdiags(ones(N1+1,1),0,N1+1,N1+1);
% add boundary conditions, 2nd order interpolation 
ddx1_ce(1,1) = 3;
ddx1_ce(1,2) = -1/3;
ddx1_ce(N1+1,N1) = -3;
ddx1_ce(N1+1,N1-1) = 1/3;
ddx1_ce = (1/dx1)*ddx1_ce(1:N1+1,1:N1);
% compute diffusive flux
diff_flux = -DM1*kron(speye(N2),ddx1_ce);
% compute advective flux
adv_flux = u1*interpx1;
flux_matrix = -(adv_flux+diff_flux);    
% projection_matrix
Pf = kron(ones(1,N2), 1/n*speye(N1+1));   
flux_matrix_bar = Pf*flux_matrix*LEC*L_bar;

aug_ddx1_ce = [zeros(N1+1,1) ddx1_ce];
aug_ddx1_ce(1,1) = 1;

aug_flux_matrix_bar = [zeros(N1+1,1) flux_matrix_bar];
aug_flux_matrix_bar(1,1) = DM1;

M = aug_flux_matrix_bar*inv(aug_ddx1_ce);

% check
disp('...checking M');
disp('2-norm')
disp(norm(full(flux_matrix_bar)-M*ddx1_ce));
disp('max diff')
disp(max(max(abs(full(flux_matrix_bar)-M*ddx1_ce))));

% eddy diffusivity 
D = M - DM1*speye(N1+1);

figure()
pcolor(log10(abs(D)))
shading flat
axis ij
colorbar
title('$D$','Interpreter','Latex')

figure()
plot(log10(abs(D(index,:))))

%% Compute Moments
D0 = zeros(N1+1,1);
D1 = zeros(N1+1,1);
D2 = zeros(N1+1,1);
for i=1:N1+1
    eta = x1_f - x1_f(i);               % eta = (y1 - x1)
    D0(i) = trapz(eta, D(i,:)/dx1);     % integrate along rows
    D1(i) = trapz(eta, eta.*D(i,:)/dx1);
    D2(i) = trapz(eta, 0.5*eta.^2.*D(i,:)/dx1);
end

% alternative way of computing moments
%D0_check = D*ones(N1+1,1);  % very little difference from D0
%D1_check = D*x1_f' - D0_check.*x1_f';
%D2_check = D*(x1_f'.*x1_f'/2) - D0_check.*(x1_f'.*x1_f'/2) - D1_check.*x1_f';

%{
figure()
hold on
box on
plot(x1_f,D0,'LineWidth',2)
plot(x1_f,D1,'LineWidth',2)
plot(x1_f,D2,'LineWidth',2)
xlim([-pi pi])
xlabel('x_1','FontSize',16)
ylabel('D^n','FontSize',16)
legend('D^0','D^1','D^2')
set(gca,'FontSize',16)
%}

% make moments symmetric
%D0 = (D0 + flipud(D0))/2;
%D1 = (D1 - flipud(D1))/2;
%D2 = (D2 + flipud(D2))/2;


%% Compute Lbar' (Lbar with molecular part subtracted off)
display('compute Lbar prime')
L_bar_p = L_bar + DM1*d2dx12_bar;

%% With forcing from manuscript
% DNS
RHS = ones(N1*N2,1); % test forcing from Mani and Park (2021)
c = L\RHS;
% compute cbar
cbar = Pc*c;

%{
% plot
figure()
contour(x1,x2,reshape(c,N1,N2)',20,'LineWidth',1.5)
colormap jet
h = colorbar;
set(gca,'FontSize',16)
set(get(h,'label'),'string','$c$','Interpreter','Latex','FontSize',20);
xlabel('$x_1$','Interpreter','Latex','FontSize',20)
ylabel('$x_2$','Interpreter','Latex','FontSize',20)

figure()
plot(x1,cbar,'o-','LineWidth',2)
xlim([-pi pi])
%ylim([0 60])
xlabel('x_1','FontSize',16)
ylabel('$\bar{c}$','Interpreter','Latex','FontSize',16)
legend('DNS')
set(gca,'FontSize',16)
%}


%% Models
% Boussinesq model 
LHS_0th = -ddx1_ec*(D0.*ddx1_ce) - DM1*d2dx12_bar; 
RHS_0th = ones(N1,1);  
cbar_0th = LHS_0th\RHS_0th;

% figure()
% hold on
% box on
% plot(x1,cbar,'k','LineWidth',2)
% plot(x1,cbar_0th,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
% xlim([-pi pi])
% ylim([0 55])
% set(gca,'FontSize',16)
% xlabel('$x_1$','Interpreter','Latex','FontSize',20)
% ylabel('$\bar{c}$','Interpreter','Latex','FontSize',20)
% leg = legend('DNS','Boussinesq model','Location','SouthEast');
% set(leg,'Interpreter','Latex')
% legend('boxoff')

%% Compare with Florian's reconstructions
directory = 'data/';
addpath(directory)
D_reconstruction = readmatrix('100N-1Nvx-2Nvy-D_2_times_21_matvecs.csv');

error()
cbar_reconstruction = Lbar_reconstruction\RHS_0th; 


figure()
hold on
box on
plot(x1,cbar,'k','LineWidth',2)
plot(x1,cbar_0th,'--','Color',[0.5 0.5 0.5],'LineWidth',2)
plot(x1,cbar_reconstruction,'LineWidth',2)
xlim([-pi pi])
%ylim([0 55])
set(gca,'FontSize',16)
xlabel('$x_1$','Interpreter','Latex','FontSize',20)
ylabel('$\bar{c}$','Interpreter','Latex','FontSize',20)
leg = legend('DNS','Boussinesq model','n191 reconstruction','Location','SouthEast');
set(leg,'Interpreter','Latex')
legend('boxoff')
