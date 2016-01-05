function [eqPt] = func_eq_pts_rolling_ball(eqNum)

%FUNC_EQ_PTS_ROLLING_BALL solves the saddle center equilibrium point for
%   the ball rolling on a surface under uniform gravity. 
%
%   [eqPt] = FUNC_EQ_PTS_ROLLING_BALL(eqNum) returns the equilibrium point
%   eqPt as (x, y) 2-tuple in the configuration space for the equilibrium
%   point indexed eqNum
%--------------------------------------------------------------------------
%   BALL ROLLING ON THE SURFACE CONVENTION:
%
%    2nd quadrant   L12(EQNUM=3)    1st quadrant
%
%
%     L23(EQNUM=4)                  L41(EQNUM=2)
%
%
%    3rd quadrant   L34(EQNUM=1)    4th quadrant
%
%--------------------------------------------------------------------------
% 
%   Shibabrat Naik (25-Nov-2015)
%   
    
    %fix the equilibrium point numbering convention here and make a
    %starting guess at the solution
    if 	eqNum == 1, 
        x0 = [0; -5];   % L34  
    elseif 	eqNum == 2, 
        x0 = [5; 0];    % L41
    elseif 	eqNum == 3, 
        x0 = [0; 5];    % L12
    elseif 	eqNum == 4, 
        x0 = [-5; 0];   % L23
    end
    
    %F(xEq) = 0 at the equilibrium point, solve using in-built function
    
    options = optimoptions('fsolve','Display','iter'); % Option to display output
    [eqPt,fval] = fsolve(@func_vec_field_eq_pt,x0,options) % Call solver

end
function F = func_vec_field_eq_pt(x)
    
    %constant parameters for the surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;

    Hx = alpha*(2*x(1)) - beta*(x(1)/(sqrt(x(1)^2 + gamma))) - xi*x(2);
    Hy = alpha*(2*x(2)) - beta*(x(2)/(sqrt(x(2)^2 + gamma))) - xi*x(1);
    
    F = [Hx/(1 + Hx^2 + Hy^2);
        Hy/(1 + Hx^2 + Hy^2)];
    
end



