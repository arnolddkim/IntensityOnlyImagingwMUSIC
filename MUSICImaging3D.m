%% MUSICImaging3D.m
%
% Compute the procedure given in Section 3 of "Intensity-only inverse
% scattering using MUSIC" by A. D. Kim and C. Tsogka. This procedure
% applies a modification of the MUSIC algorithm for the B data matrix
% computed by ForwardProblem3D.m.
%
% Written by C. Tsogka on 4/17/2019. Revised by A. D. Kim on 5/2/2019.

clear;

%% Load the forward data

load ForwardData;

%% STEP 1: Compute the singular value decomposition of the B matrix

% extract rows 2 to the end of the B matrix

Bdata = B(2:end,:);
 
% compute the SVD of B

[ U, Sigma, V ] = svd( Bdata );

% find the significant singular values

sig     = diag( Sigma );
nsignal = length( find( sig / max( sig ) > 1.e-8 ) );

%% STEP 2: Compute the projection

Projection = eye( N_Tx - 1 ) - U(:,1:nsignal) * U(:,1:nsignal)';

%% STPE 3: Compute eta_k = || P a_k ||

% compute the imaging region

[ X, Y, Z ] = meshgrid( -10 * lambda : 0.25 * lambda : 10 * lambda );

% compute meshgrid of indices

Ngrid = length( X );

[ indx, jndx ] = ndgrid( ( 2 : N_Tx ), ( 1 : Ngrid^3 ) );

% compute the projection of the A matrix of incident fields

sx = dir(:,1);
sy = dir(:,2);
sz = dir(:,3);

% compute P * a_k

Pa_k = Projection * exp( 1i * 2 * pi / lambda * ( sx(indx) .* X(jndx) ...
    + sy(indx) .* Y(jndx) + sz(indx) .* Z(jndx) ) );

% compute eta

eta = sqrt( sum( abs( Pa_k ).^2, 1 ) );

%% STEP 4: Compute I_k = min(eta_k)/eta_k

Image_MUSIC = reshape( min(eta) ./ eta, Ngrid, Ngrid, Ngrid );

%% set the figure parameters

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% plot isosurfaces of the image

figure(1)

p = patch( isosurface( X / lambda, Y / lambda, Z / lambda, Image_MUSIC, 0.1 ) ); 
isonormals( X / lambda, Y / lambda, Z / lambda, Image_MUSIC, p );
set(gca,'Zdir','reverse')
p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = '0.4';

%% plot the true scattering objects

hold on

for nobjs = 1:Nobjects

     x = rbdy(:,1,nobjs) / lambda;
     y = rbdy(:,2,nobjs) / lambda;
     z = rbdy(:,3,nobjs) / lambda;

    [ center, radii, evecs, v2, chi2 ] = ellipsoid_fit( [ x y z ], '' );

    mind = min( [ x y z ] );
    maxd = max( [ x y z ] );
    nsteps = 50;
    step = ( maxd - mind ) / nsteps;
    
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), ...
        linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), ...
        linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );

    Ellipsoid = v2(1) *x.*x +   v2(2) * y.*y + v2(3) * z.*z + ...
              2*v2(4) *x.*y + 2*v2(5)*x.*z + 2*v2(6) * y.*z + ...
              2*v2(7) *x    + 2*v2(8)*y    + 2*v2(9) * z;

    p2 = patch( isosurface( x, y, z, Ellipsoid, -v2(10) ) );

    set(gca,'Zdir','reverse')
    p2.FaceColor = 'blue';
    p2.EdgeColor = 'none';
    p2.FaceAlpha = '0.8';

end

hold off
 
daspect([1 1 1])
view(3); 
axis([-10 10 -10 10 -10 10]);
camlight %(-80,-10) 
lighting gouraud
grid on
grid minor
    
xlabel( '$x / \lambda$', 'Interpreter', 'LaTeX' );
ylabel( '$y / \lambda$', 'Interpreter', 'LaTeX' );
zlabel( '$z / \lambda$', 'Interpreter', 'LaTeX' );

%for normal view
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,-1,1],'Rotation',-35)

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,1,1],'Rotation',20)