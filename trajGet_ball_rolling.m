function [x,t] = trajGet_ball_rolling(x0,tb,tf,param,OPTIONS)

%        [x,t] = trajGet_ball_rolling(x0,tb,tf,param,OPTIONS);
%
% This is an integrator for a point in phase space
%
%         d x(t)
%        ------ =  f(x)
%          dt
%
%
% Shane Ross (revised 10.13.03)
% Shibabrat Naik (modified for ball rolling on surface: 2015-Dec-05)

MODEL = 'ball_rolling2dof';

if nargin==4;
   OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12);  % lower accuracy
end

TSPANtb = [-1.e-12 tb ] ;
TSPANtf = [ 0      tf ] ;

x=[];

if tb==0,
  [t,x]     = ode113(MODEL,TSPANtf,x0,OPTIONS);
elseif tf==0,
  [t,x]     = ode113(MODEL,TSPANtb,x0,OPTIONS);
else,
  [tt1,xx1] = ode113(MODEL,TSPANtb,x0,OPTIONS);
  [tt2,xx2] = ode113(MODEL,TSPANtf,x0,OPTIONS);

  x=[flipud(xx1);xx2];
  t=[flipud(tt1);tt2];
end

end