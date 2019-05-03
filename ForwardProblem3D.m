%% ForwardProblem3D.m
%
% Compute scattering by scattering objects due to several incident plane
% waves at different angles of incidence using the method of fundamental 
% solutions (MFS).
%
% The scattered fields on a measurement plane are collected and coverted
% into the interferometric data used for imaging in the paper entitled,
% "Intensity-only inverse scattering using MUSIC" by A. D. Kim and C.
% Tsogka (2019).
%
% This code requires the following other codes.
%
% * Ellipsoid.m
% * ComputeMFSCoeffs.m
% * ComputeMeasurements.m
% 
% The key result is the B matrix given in Eq. (6) of Kim and Tsogka.
%
% Written by A. D. Kim on 2/16/2019 and modified most recently on 4/17/2019

clear;

%% PHYSICAL PARAMETERS

% wavelength in nanometers

lambda = 632;

% refractive indices

n0 = 1.00; % exterior
n1 = 1.40; % interior

% wavenumbers 

k0 = 2 * pi * n0 ./ lambda; % exterior
k1 = 2 * pi * n1 ./ lambda; % interior

%% BOUNDARY/INTERIOR/EXTERIOR POINTS

N = 512; % number of points

% set the number of scattering objects

Nobjects = 3;

% set initial locations for scattering objects

R0 = 8 * [ 0 0 0; 0 1 0; -1 0 0 ];
    
% randomly perturb those locations

R0 = lambda * ( R0 + 0.1 * randn( Nobjects, 3 ) );

% allocate memory for the MFS points for each object

rbdy = zeros( N, 3, Nobjects );
rint = zeros( N, 3, Nobjects );
rext = zeros( N, 3, Nobjects );
nu   = zeros( N, 3, Nobjects );

% compute the locations of the scattering objects

for nobjs = 1 : Nobjects
   
    % compute the reference object
    
    a = lambda;
    b = 1.50 * a;
    c = 0.75 * a;
    
    [ nu0, rbdy0, rint0, rext0 ] = Ellipsoid( a, b, c, N );
    
    % translate the object by R0
    
    rbdy0 = R0(nobjs,:) + rbdy0;
    rint0 = R0(nobjs,:) + rint0;
    rext0 = R0(nobjs,:) + rext0;
    
    % rotate the object
    
    theta0 = pi * rand(1);
    phi0   = 2 * pi * rand(1);
    
    Rot = [ cos(theta0)*cos(phi0), -sin(phi0), sin(theta0)*cos(phi0);
            cos(theta0)*sin(phi0),  cos(phi0), sin(theta0)*sin(phi0);
           -sin(theta0),            0,         cos(theta0) ];

    rbdy(:,:,nobjs) = rbdy0 * Rot';
    rint(:,:,nobjs) = rint0 * Rot';
    rext(:,:,nobjs) = rext0 * Rot';
    nu(:,:,nobjs)   = nu0   * Rot';
    
    figure(100);
    plot3( rbdy(:,1,nobjs) / lambda, rbdy(:,2,nobjs) / lambda, ...
        rbdy(:,3,nobjs) / lambda, '.' );
    xlabel( '$x / \lambda$', 'Interpreter', 'LaTeX' );
    ylabel( '$y / \lambda$', 'Interpreter', 'LaTeX' );
    zlabel( '$z / \lambda$', 'Interpreter', 'LaTeX' );
    title( 'Location of ellipsoids', 'Interpreter', 'LaTeX' );
    hold on;
        
end

grid;
grid minor;
hold off;
  
%% SOURCES: set the angles of incidence for the illuminating plane waves

Mtheta = 25;
Mphi   = 25;

% define polar angles

theta = linspace( 0, pi / 2, Mtheta + 2 );
theta = theta(2:end-1);

% define azimuthal angles

phi   = linspace( 0, 2 * pi, Mphi );

% define direction vectors

[ THETA, PHI ] = meshgrid( theta, phi );
THETA = THETA(:);
PHI   = PHI(:);
dir   = [ 0 0 1; sin(THETA).*cos(PHI) sin(THETA).*sin(PHI) cos(THETA) ];
N_Tx  = length( dir );

%% RECEIVERS: set the array of detectors

Rx_mesh = ( -200 * lambda : 10 * lambda : 200 * lambda );
N_Rx    = length( Rx_mesh );

[ X_Rx, Y_Rx ] = meshgrid( Rx_mesh );
X_Rx = X_Rx(:);
Y_Rx = Y_Rx(:);
Z_Rx = 200 * lambda;

%% SCATTERED FIELD MEASUREMENTS

% allocate memory for the scattered field data

Udata = zeros( N_Tx, N_Rx^2 );

% loop over each source and collect the scattered field data (uses Matlab
% Parallel Toolbox to speed up this computation -- if that toolbox is not
% available, change "parfor" to "for" in the line below)

parfor n = 1 : N_Tx

    % compute the MFS coefficients

    [ ~, csca ] = ComputeMFSCoeffs( rbdy, rint, rext, nu, dir(n,:), ...
        k0, k1, N, Nobjects );

    % compute the scattered field measurements
    
    Udata(n,:) = ComputeMeasurements( k0, csca, rext, X_Rx, Y_Rx, Z_Rx, ...
        Nobjects );
    
    % inform the user on progress    
    
    disp( [ '  Computed source ', num2str(n), ' out of ', num2str(N_Tx) ] );
    
end

%% FORM INTERFEROMETRIC DATA MATRIX B

% allocate memory for the B matrix

B = zeros( N_Tx, N_Rx^2 );

% set the index for the reference field (normal incidence)

istar = 1;

for nr = 1 : N_Rx^2
   
    B(:,nr) = conj( Udata(istar,nr) ) * Udata(:,nr);
    
end

%% SAVE THE DATA

save ForwardData B N_Rx N_Tx Nobjects rbdy dir lambda;