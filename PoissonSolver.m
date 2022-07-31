%% Contents
% We present the new Poisson solver.

%% Simulation Parameters (grid, etc.)
clearvars
format longg

Lx = 30; Ly = 30; Lz = 60; % Spatial extent of box in Fermi
dx = 1; dy = dx; dz = dx; % Lattice spacing in Fermi
Nx = Lx/dx; Ny = Ly/dy; Nz = Lz/dz;  dV = dx*dy*dz; 
N = max([Nx,Ny,Nz]); L = max([Lx,Ly,Lz]); N3 = 3*N; 
xaxis = dx*(0:Nx-1); yaxis = dy*(0:Ny-1); zaxis = dz*(0:Nz-1); % Spatial grid
[X,Y,Z] = ndgrid(xaxis, yaxis, zaxis);
e2=197.3269631 / 137.035999679; % Squared charge of proton

% Allowed array sizes for new solver
sizes = [2, 4, 6, 8, 10, 12, 16, 18, 20, 24, 30, 32, 36, 40, 48, 50, 54, 60, 64, 72, 80, 90, 96, 100, 108, 120, 128, 144, 150, 160, 162, 180, 192, 200, 216, 240, 250, 256, 270, 288, 300];

solver = 2; % 1: old solver. 2: new solver.
%% Define gaussian test charge distribution

sig = 1.87; % Standard deviation of charge distribution
numimg = 2; % Number of image charges

rho = zeros(Nx,Ny,Nz);
% Add image charges on all sides
for imagex = -numimg:numimg
    for imagey = -numimg:numimg
        for imagez = -numimg:numimg
            rho = rho + exp(-((X-Lx/2+imagex*Lx).^2+...
                              (Y-Ly/2+imagey*Ly).^2+...
                              (Z-Lz/2+imagez*Lz).^2)/(2*sig^2));
        end
    end
end
rho = rho/((2*pi*sig^2)^(3/2)); % Normalize 
disp("Total charge:")
disp(sum(rho,"all")*dV)

%% Solve the Poisson equation

tic

% Old solver
if solver == 1 
    
    % Radius of Coulomb kernel.
    D = sqrt(3)*L;
    
    % Embed rho in larger box
    rho3 = zeros(N3,N3,N3);
    rho3(1:Nx,1:Ny,1:Nz) = rho;
    
    % Define momenta and Coulomb kernel in larger box
    k2 = make_ks(N3,N3,N3,dx,dy,dz);
    coulomb_kernel = get_coul_ker(k2,D,e2);
    
    % Compute potential in larger box
    phi3 = ifftn(coulomb_kernel.*fftn(rho3));
    phi = real(phi3(1:Nx,1:Ny,1:Nz));
    
% New solver
elseif solver == 2 
    
    % Instead of extending to "N3" grid points, we try to choose smallest
    % possible grid. Then we increase array dimensions until they have 
    % small prime factors for efficient FFT.
    
    % Radius of Coulomb kernel
    D = sqrt(Lx^2+Ly^2+Lz^2);
    
    % Initial guess for grid extension size
    dNx0 = 2*ceil(ceil(D/dx)/2);
    dNy0 = 2*ceil(ceil(D/dy)/2);
    dNz0 = 2*ceil(ceil(D/dz)/2);

    % Choose N_+dN_ to have small primes in prime factorization
    % to increase FFT efficiency. Factorizations are 
    % precomputed up to size 300, stored in "sizes."
    if Nx+dNx0<=300
        dNx = sizes(find(Nx+dNx0<=sizes,1))-Nx;
    else
        disp("Extended axis larger than 300 points")
        dNx = dNx0;
    end
    if Ny+dNy0<=300
        dNy = sizes(find(Ny+dNy0<=sizes,1))-Ny;
    else
        disp("Extended axis larger than 300 points")
        dNy = dNy0;
    end
    if Nz+dNz0<=300
        dNz = sizes(find(Nz+dNz0<=sizes,1))-Nz;
    else
        disp("Extended axis larger than 300 points")
        dNz = dNz0;
    end

    % Embed rho in larger box
    rhonew = zeros(Nx+dNx,Ny+dNy,Nz+dNz);
    rhonew(1:Nx,1:Ny,1:Nz) = rho;
    
    % Define momenta and Coulomb kernel in larger box
    k2 = make_ks(Nx+dNx,Ny+dNy,Nz+dNz,dx,dy,dz);
    coulomb_kernel = get_coul_ker(k2,D,e2);
    
    % Compute potential in larger box
    phinew = ifftn(coulomb_kernel.*fftn(rhonew));
    phi = phinew(1:Nx,1:Ny,1:Nz);    
end

disp("Solver "+solver+": elapsed time " + toc)

disp("Comparing Coulomb energy with analytical value for gaussian test case: dE/E = ")
gauss_coul_en = e2/(2*sig*sqrt(pi));
disp( abs( (1/2)*sum(rho.*phi,"all")*dV - gauss_coul_en )/gauss_coul_en )
%% Helper functions

% Define momentum lattice, return k^2
function k2 = make_ks(Nx,Ny,Nz,dx,dy,dz)
    kx = 2*pi/(Nx*dx)*[0:Nx/2-1,-Nx/2:-1];
    ky = 2*pi/(Ny*dy)*[0:Ny/2-1,-Ny/2:-1];
    kz = 2*pi/(Nz*dz)*[0:Nz/2-1,-Nz/2:-1];
    [kxgrid,kygrid,kzgrid] = ndgrid(kx,ky,kz);
    k2 = kxgrid.^2+kygrid.^2+kzgrid.^2;
end

% Define coulomb kernel
function coul_ker = get_coul_ker(k2,D,e2)
    coul_ker = 4*pi*e2*(1-cos(D*sqrt(k2)))./k2;
    coul_ker(1,1,1) = 2*e2*pi*D^2;
end

