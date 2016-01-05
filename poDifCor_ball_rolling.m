function [x0po,t1] = poDifCor_ball_rolling(x0)

%        [x0po,t1] = poDifCor_ball_rolling(x0) ;
%
% This is the differential correction routine to create a periodic
% orbit (po) about an equilibrium point. It keeps the initial  x
% value constant and varies the y value.
%
% output: xOpo  = initial state on the po (on the negative x-axis) 
%           t1  = half-period of po
%
% input:  x0    = first guess of initial state on the po 
%
%-----------------------------------------------------------------------
%
% Shane Ross (revised 2009-July-07)
% Shibabrat Naik (modified for ball rolling on the surface: 2015-Dec-20)

% tolerances for integration and perpendicular crossing of x-axis

% MAXdxdot1 = 1.e-8 ; RelTol = 3.e-10; AbsTol = 1.e-10; 
MAXdxdot1 = 1.e-12 ; RelTol = 3.e-14; AbsTol = 1.e-14; 

MAXattempt = 25;     	% maximum number of attempts before error is declared

dxdot1 	   = 1;         % to start while loop
dydot1 	   = 1;         % to start while loop
attempt    = 0;         % begin counting number of attempts
% y0(attempt) = 1;

%   constant parameters for the rolling surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
%     H0 = 12.065;
    g = 981;
    
    Hx = @(x,y)(2*alpha*x - beta*(x/(sqrt(x^2 + gamma))) - xi*y);
    Hy = @(x,y)(2*alpha*y - beta*(y/(sqrt(y^2 + gamma))) - xi*x);
    Hxx = @(x,y)(2*alpha - (beta*gamma)/((x^2 + gamma)^(3/2)) );
    Hyy = @(x,y)(2*alpha - (beta*gamma)/((y^2 + gamma)^(3/2)) );
    Hxy = @(x,y)(-xi);
    
    h = @(x,y,vx,vy)( (5*g)/7 + Hxx(x,y)*vx^2 + ...
        2*Hxy(x,y)*vx*vy + Hyy(x,y)*vy^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2);
       

while abs(dydot1) > MAXdxdot1
    
	if attempt > MAXattempt
		ERROR = 'Maximum iterations exceeded' ;
		disp(ERROR) ;
		break
    end
    
    y0 = x0;
    % Find first half-period crossing event
    MODEL = 'ball_rolling2dof';  
    TSPAN = [0 10] ;        % allow sufficient time for the half-period crossing event         
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events','on'); 
    [tt,xx,t1,xx1,i1] = ode113(MODEL,TSPAN,x0,OPTIONS);

	x1 = xx1(end,1); 
    y1 = xx1(end,2); 
	dxdot1 = xx1(end,3); 
    dydot1  = xx1(end,4); 
%     plot3(xx(:,1),xx(:,2),xx(:,3),'-r');hold on
%     plot3(x1,y1,dxdot1,'xb');
   
    
    % Compute the state transition matrix from the initial state to
	% the final state at the half-period event crossing
      
    % Events option not necessary anymore
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
	[x,t,phi_t1,PHI] = stateTransitionMatrix_ball_rolling(x0,t1,OPTIONS) ;

	attempt=attempt+1 ;
	ATTEMPT = sprintf('::poDifCor : iteration %d',attempt) ;
	disp(ATTEMPT) ;
    
    %	show = 1 to plot successive orbit  (default=0)
    show = 1 ; % set to 1 to plot successive approximations
    if show==1,
        plot(x(:,2),x(:,3),'.-'); 
        hold on;
        m = length(x) ;
        plot(x(1,2),x(1,3),'b*');
        plot(x(m,2),x(m,3),'bo');
    %         axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);
    %         axis equal
        xlabel('$y(t)$','interpreter','latex','fontsize',24);
        ylabel('$v_x(t)$','interpreter','latex','fontsize',24);
        set(gca,'fontsize',18)
        pause(0.01) ;
    end


%=========================================================================
% differential correction and convergence test, adjust according to
% the particular problem

    %compute acceleration values for correction term
%     vxdot1 = -x1 + 2*x1*y1;
%     vydot1 = -R^2*y1 + 0.5*R^2*x1^2;
    vxdot1 = -h(x1,y1,dxdot1,dydot1)*Hx(x1,y1);
    vydot1 = -h(x1,y1,dxdot1,dydot1)*Hy(x1,y1);
    
    %correction to the initial y0
%     y0(attempt) = 1/(phi_t1(3,2) - phi_t1(4,2)*vxdot1*(1/vydot1))*dxdot1;
	y0(attempt) = 1/(phi_t1(4,2) - phi_t1(3,2)*vydot1*(1/vxdot1))*dydot1;
	
    x0(2) = x0(2) - y0(attempt);
    
end

x0po=x0;
t1 = t1(end);

end


