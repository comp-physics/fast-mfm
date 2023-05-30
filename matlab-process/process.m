close all; clear all; clc;
%%% Solver for the inhomogeneous model problem in Mani and Park (2021) %%%
% staggered mesh using 2nd order finite difference
% scalar c is stored at cells
% velocities u1 and u2 are stored at faces


directory = 'data/';
addpath(directory)

directory = '/Users/spencer/Downloads/multivortex/case-2000N-3Nvx-2Nvy/';
addpath(directory)

directory = '/Users/spencer/research/docs/paper_mfm/data/laminar_channel/matrices/';
% directory = '/Users/spencer/Downloads/bigmat/'
addpath(directory)

%% Parameters
DM1 = 0.05;
DM2 = 1;

% number of mesh points
N1 = 2000; %40;   
N2 = N1/2; %10;    

Nvx = 3;
Nvy = 2;

if mod(Nvx,2) == 0
    phaseskipX = pi;
else
    phaseskipX = 0;
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
L = ddx1*(u1*interpx1) + ddx2*(u2*interpx2) - DM1*d2dx12 - DM2*d2dx22; 

%% P Matrix   
% averaged in x2 direction
n = N2;
Pc = kron(ones(1,N2), 1/n*speye(N1));

%% E Matrix
Ec = n*Pc';

%% Eddy-Diffusivity 
% ddx1 cell center to edge
ddx1_ce = spdiags(-ones(N1+1,1),-1,N1+1,N1+1)+spdiags(ones(N1+1,1),0,N1+1,N1+1);
% add boundary conditions, 2nd order interpolation 
ddx1_ce(1,1) = 3;
ddx1_ce(1,2) = -1/3;
ddx1_ce(N1+1,N1) = -3;
ddx1_ce(N1+1,N1-1) = 1/3;
ddx1_ce = (1/dx1)*ddx1_ce(1:N1+1,1:N1);


RHS_0th = ones(N1,1);  

D = readmatrix('2000N1-1000N2-3Nvx-2Nvy-D.csv');
trueD = transpose(interp1)*D*interp1;
trueL = -ddx1_ec*(D*ddx1_ce) - DM1*d2dx12_bar;
cbar = trueL\RHS_0th; 

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

%% Models
% Boussinesq model 
LHS_0th = -ddx1_ec*(D0.*ddx1_ce) - DM1*d2dx12_bar; 
cbar_0th = LHS_0th\RHS_0th;

%% Compare with Florian's reconstructions
% matrices = { ...
%     '100N-1Nvx-2Nvy-D_2_times_21_matvecs.csv', ... 
%     '100N-1Nvx-2Nvy-D_2_times_29_matvecs.csv', ...
%     '100N-1Nvx-2Nvy-D_2_times_26_matvecs.csv', ...
%     '100N-1Nvx-2Nvy-D_2_times_33_matvecs.csv' ...
%     };

matrices = {...
%     '2000N-3Nvx-2Nvy-D_2_times_11_matvecs.csv'
%     '2000N-3Nvx-2Nvy-D_2_times_15_matvecs.csv' ...
      '2000N-3Nvx-2Nvy-D_rand_lr_n_meas_26.csv' ...
%     '2000N-3Nvx-2Nvy-D_2_times_21_matvecs.csv' ...
%  '2000N-3Nvx-2Nvy-D_2_times_61_matvecs.csv', ...
%   '2000N-3Nvx-2Nvy-D_2_times_81_matvecs.csv' 
};

% matrices = {...
%      '2000N-1Nvx-1Nvy-D_2_times_11_matvecs.csv' ...
%      '2000N-1Nvx-1Nvy-D_2_times_21_matvecs.csv' ...
%    '2000N-1Nvx-1Nvy-D_2_times_21_matvecs.csv' ...
%    '2000N-1Nvx-1Nvy-D_2_times_41_matvecs.csv' ...
%    '2000N-1Nvx-1Nvy-D_2_times_61_matvecs.csv', ...
%    '2000N-1Nvx-1Nvy-D_2_times_81_matvecs.csv' 
% };

for ii = 1:length(matrices)
    disp(matrices{ii})
    D_reconstruction = readmatrix(matrices{ii});
    disp('Dreconstruct error')
    norm(D_reconstruction - D)/norm(D)

    newD = transpose(interp1)*D_reconstruction*interp1;
    newL = -ddx1_ec*(D_reconstruction*ddx1_ce) - DM1*d2dx12_bar;
    cbar_reconstruction = newL\RHS_0th; 

    disp('cbar err')
    disp(norm(cbar_reconstruction-cbar)/norm(cbar))
end


csvwrite('data/2000N-3Ny-2Nx-exact.csv',[x1;cbar']');
csvwrite('data/2000N-3Ny-2Nx-bous.csv',[x1;cbar_0th']');
csvwrite('data/2000N-3Ny-2Nx-rand-26.csv',[x1;cbar_reconstruction']');

figure()
hold on
box on
plot(x1,cbar,'k','LineWidth',2.5)
plot(x1,cbar_0th,'--','Color',[0.72 0.26 0.06],'LineWidth',2.5)
% plot(x1,cbar_reconstruction,'Color',[0.18 0.42 0.41],'LineWidth',2)
plot(x1(1:100:2000),cbar_reconstruction(1:100:2000),'-s','MarkerSize',12,...
        'MarkerFaceColor',[0.18 0.42 0.41], ...
    'MarkerEdgeColor','black', ...
    'LineStyle', 'none' )

xlim([-pi pi])
set(gca,'FontSize',22)
xlabel('$x$','Interpreter','Latex','FontSize',22)
ylabel('$\bar{c}$','Interpreter','Latex','FontSize',22)
leg = legend('DNS','Boussinesq model','N=21 reconstruction (ours)','Location','SouthEast');
set(leg,'Interpreter','Latex')
legend('boxoff')



% figure()
% pcolor(log10(abs(newL)))
% shading flat
% axis ij
% colorbar
% title('Lnew','Interpreter','Latex')

% figure()
% pcolor(log10(abs(D_reconstruction)))
% shading flat
% axis ij
% colorbar
% title('$Dreconstruct$','Interpreter','Latex')

% error()
