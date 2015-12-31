function test_periodic_orbit

%TEST_PERIODIC_ORBIT functions integrates the periodic orbit initial
%condition

    y0= [-0.332463579671972; ...
        6.122156347632109; ...
                        0; ...
                        0];
    
                    
    tspan = [0 1];
    RelTol = 1e-12;
    AbsTol = 1e-10;
    options = odeset('RelTol',1e-12,'AbsTol',1e-10,...
                 'OutputFcn',@odephas2,'Events',@events);

    MAXattempt = 25;     	% maximum number of attempts before error is declared

    xdot1 	   = 1;   	% to start while loop
    attempt    = 0;		% begin counting number of attempts
    MAXdxdot1 = 1.e-10 ;
    
    attempt = 0;
    while abs(xdot1) > MAXdxdot1
	
    if attempt > MAXattempt
		ERROR = 'Maximum iterations exceeded' ;
		disp(ERROR) ;
		break
	end

    [t,y,te,ye,ie] = ode113(@f,tspan,y0,options)

    plot3(y(:,1),y(:,2),y(:,3),'-b');%,ye(:,1),ye(:,2),'o');
    ylabel ('y(t)')
    xlabel ('x(t)')
    grid on 
    
    x1     = ye(end,1) ; 
    y1     = ye(end,2) ; 
	xdot1 =  ye(end,3) ; 
    ydot1  = ye(end,4) ; 
    t1 = 2*te(end);

    
    % Compute the state transition matrix from the initial state to
	% the final state 
    % Events option not necessary anymore
    OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
	[x,t,phi_t1,PHI] = stateTransitionMatrix_ball_rolling(y0,t1,1,OPTIONS) ;

	attempt=attempt+1 ;
	ATTEMPT = sprintf('::poDifCor : iteration %d',attempt) ;
	disp(ATTEMPT) ;
    
    %=========================================================================
    % differential correction and convergence test, adjust according to
    % the particular problem

    %compute acceleration values for correction term
    vxdot1 = -h(x1,y1,xdot1,ydot1)*Hx(x1,y1);
    vydot1 = -h(x1,y1,xdot1,ydot1)*Hy(x1,y1);
    
    %correction to the initial state, y component 
    yCorr(attempt) = 1/(phi_t1(3,2) - phi_t1(4,2)*vxdot1*(1/vydot1))*xdot1;
    
	y0(2) = y0(2) - yCorr(attempt)
    
    end

% --------------------------------------------------------------
function dydt = f(t,y)
    
    %constant parameters for the rolling surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
    H0 = 12.065;
    g = 918;

    %function for the rolling surface 
%     H = @(x,y)(alpha*(x^2 + y^2) - ...
%         beta*(sqrt(x^2 + gamma) + sqrt(y^2 + gamma)) - ...
%         xi*x*y + H0); 
    Hx = @(x,y)(2*alpha*x - beta*(x/(sqrt(x^2 + gamma))) - xi*y);
    Hy = @(x,y)(2*alpha*y - beta*(y/(sqrt(y^2 + gamma))) - xi*x);
    Hxx = @(x,y)(2*alpha - (beta*gamma)/((x^2 + gamma)^(3/2)) );
    Hyy = @(x,y)(2*alpha - (beta*gamma)/((y^2 + gamma)^(3/2)) );
    Hxy = @(x,y)(-xi);
    
    h = @(x,y,vx,vy)( (5*g)/7 + Hxx(x,y)*vx^2 + ...
        2*Hxy(x,y)*vx*vy + Hyy(x,y)*vy^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2);
    
    dydt = zeros(4,1);
        
    dydt(1) = y(3);
    dydt(2) = y(4);
    dydt(3) = -h(y(1),y(2),y(3),y(4))*Hx(y(1),y(2));
    dydt(4) = -h(y(1),y(2),y(3),y(4))*Hy(y(1),y(2));
    
end   % End nested function f
% --------------------------------------------------------------
function [value,isterminal,direction] = events(t,yt)
% Locate the time when the object returns closest to the
% initial point y0 and starts to move away; stop integration.
% Also locate the time when the object is farthest from the 
% initial point y0 and starts to move closer.
% 
% The current distance of the body is
% 
%   DSQ = (y(1)-y0(1))^2 + (y(2)-y0(2))^2 
%       = <y(1:2)-y0(1:2),y(1:2)-y0(1:2)>
% 
% A local minimum of DSQ occurs when d/dt DSQ crosses zero 
% heading in the positive direction. Compute d(DSQ)/dt as
% 
%  d(DSQ)/dt = 2*(y(1:2)-y0(1:2))'*dy(1:2)/dt = ... 
%                 2*(y(1:2)-y0(1:2))'*y(3:4)
% 

if abs(t) > 1e-2
    isterminal = 1;  %terminate after waiting for a short time
else
    isterminal = 0;
end

% dDSQdt = 2 * ((yt(1:2)-y0(1:2))' * yt(3:4));
% value = [dDSQdt; dDSQdt];
% isterminal = [1; 0];            % Stop at local minimum
% direction = [1; -1];            % [local minimum, local maximum]
value = yt(4);
direction = 0;

end   % End nested function events


end