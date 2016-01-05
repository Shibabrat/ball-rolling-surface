function [stateVarsDot,termEvent, direcEvent] = ...
    ball_rolling2dof(t, stateVars, flag)

%BALL_ROLLING2DOF defines the odes for the ball rolling on a surface
%   problem and has event option for one of the phase space
%   coordinatesusing a flag.
%
%   [stateVarsDot,termEvent, direcEvent] = ball_rolling2dof(t, stateVars, flag)
%

    %constant parameters for the rolling surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
    H0 = 12.065;
    g = 981;

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
       
if (nargin < 3 || isempty(flag)) %without event 

    stateVarsDot = zeros(4,1);
    x = stateVars(1); 
    y = stateVars(2);
    vx = stateVars(3);
    vy = stateVars(4);
        
    stateVarsDot(1) = vx;
    stateVarsDot(2) = vy;
    stateVarsDot(3) = -h(x,y,vx,vy)*Hx(x,y);
    stateVarsDot(4) = -h(x,y,vx,vy)*Hy(x,y);
    
else
   
    switch lower(flag)      %with event
        case 'events'
            if abs(t) > 1e-2
                isterminal = 1;  %terminate after waiting for a short time
            else
                isterminal = 0;
            end
            
            direction = 0; %0: all directions of crossing
            
            stateVarsDot = stateVars(3);
    
            termEvent = isterminal;
            direcEvent = direction;            
    end
    
end

end