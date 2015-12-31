function [x0poGuess,TGuess] = poGuessLinear_ball_rolling(eqNum,Ax)

%        [x0poGuess,TGuess] = poGuessLinear_ball_rolling(param,eqNum,Ax) ;
%
% Uses a small displacement from the equilibrium point (in a direction 
% on the collinear point's center manifold) as a first guess for a planar 
% periodic orbit (called a Lyapunov orbit in th rest. three-body problem).
%
% The initial condition and period are to be used as a first guess for
% a differential correction routine.
%
% output:
% x0poGuess = initial state on the periodic orbit (first guess)
%           = [ x 0  0 yvel]  , (i.e., perp. to x-axis and in the plane)
% TGuess    = period of periodic orbit (first guess)
%
% input:
% eqNum = the number of the equilibrium point of interest
% Ax    = nondim. amplitude of periodic orbit (<< 1) 
%
%----------------------------------------------------------------------------
%
% Shane Ross (revised 2.17.04)
% Shibabrat Naik (modified for ball rolling problem: 2015-Dec-26)

% eqPos = eqPointLoc(G,eqNum);  % position space location of equil. point
eqPos = func_eq_pts_rolling_ball(eqNum);
ep = [eqPos' 0 0];              % phase space location of equil. point

% Get the eigenvalues and eigenvectors of Jacobian of ODEs at equil. point
[Es,Eu,Ec,Vs,Vu,Vc]=eqPointEig_ball_rolling(ep);

l = abs(imag(Ec(1))) ;  

Df = func_get_Jacobian(ep);
k2 = -(l^2 + Df(3,1))/(Df(3,2));

x0poGuess       = zeros(4,1) ;
x0poGuess(1)	= ep(1) - Ax ;
x0poGuess(2)	= ep(2) - Ax*k2 ;

TGuess = 2*pi/l ;

end
