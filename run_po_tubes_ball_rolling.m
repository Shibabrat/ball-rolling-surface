%   SCRIPT to compute periodic orbits, tube manifolds for the ball rolling
%   on a surface, 2DOF problem
%==========================================================================
%   Note: It closely follows the dynamics_CR3BP code by Shane Ross and the
%   theory in chapter 4 of Dynamical Systems, the Three-Body Problem and
%   Space Mission Design W.S. Koon, M.W. Lo, J.E. Marsden and S.D. Ross
%   ISBN 978-0-615-24095-4 Marsden Books, 2008.
%
%   http://www.shaneross.com/books
%==========================================================================

%   Setting up parameters
N = 4;          %dimension of phase space
g = 981;        %Acceleration due to gravity
eqNum = 2;      %3: corresponds to the saddle eq pt between 1st and 2nd quadrant    

%   first two amplitudes for continuation procedure to get p.o. family
nFam = 300;
Ax1  = 2.e-4; % initial amplitude (1 of 2)
Ax2  = 2*Ax1; % initial amplitude (2 of 2)

tic;

%   get the initial conditions and periods for a family of periodic orbits
[x0po,TPOFam] = poFamGet_ball_rolling(eqNum,Ax1,Ax2,nFam) ; 

poFamRuntime = toc

x0podata = [x0po, TPOFam];

%%
eSaddle = 3304.24;      %energy of the saddle eq pt, for sake of reference
eTarget = 3504.24; 

%   target specific periodic orbit 
fileName = 'x0po_T_energy_case1_L41.txt';
PODATAFILE = sprintf('Loading the periodic orbit family from data file %s',fileName); 
disp(PODATAFILE)
x0podata = importdata(fileName);

[x0poTarget,TTarget] = poBracketEnergy_ball_rolling(eTarget,x0podata); 

%%
% get and plot unstable manifold (negative branch) of a p.o. specifically,
% plot 40 (roughly) equally spaced trajectories on this `tube' manifold
% whose initial conditions are displaced 1.e-6 from the p.o. in phase space
n_mfd_traj = 100;
x0po= importdata('./data/saddle-L12/x0po_T_energy_case2.txt');
TPOFam = x0po(:,5);
ePOFam = x0po(:,end);
nMed = size(x0po,1);

tic;

% [xW,x0W] = poManiLocal_ball_rolling(x0po(nMed,1:4),TPOFam(nMed),0,1,1, ...
%     1e-6,1,4*TPOFam(nMed),n_mfd_traj);
[xW,x0W] = poManiLocal_ball_rolling(x0po(nMed,1:4),TPOFam(nMed),0,-1,1, ...
    1e-6,1,4*TPOFam(nMed),n_mfd_traj);

maniRuntime = toc

title(['Energy of trajectories on the tube:',num2str(ePOFam(nMed))])

energyTube = ePOFam(nMed)








