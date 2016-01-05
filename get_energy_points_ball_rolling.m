function e = get_energy_points_ball_rolling(orbit)

%GET_ENERGY_POINTS_BALL_ROLLING computes the energy of an input orbit
%   (represented as M x N with M time steps and N = 4, dimension of phase
%   space for the model) for a ball rolling on a surface under uniform
%   gravity.
%
%   e = GET_ENERGY_POINTS_BALL_ROLLING(orbit) returns energy e of the orbit
%   in the phase space. orbit can be different initial conditions for the
%   periodic orbit of different energy.
%
%   Shibabrat Naik (Revised: 05-Jan-2016)
%

    x = orbit(:,1);
    y = orbit(:,2);
    vx = orbit(:,3);
    vy = orbit(:,4);
    
    %constant parameters for the rolling surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
    H0 = 12.065;
    g = 981;

    %function for the rolling surface 
    H = @(x,y)(alpha*(x.^2 + y.^2) - ...
        beta*(sqrt(x.^2 + gamma) + sqrt(y.^2 + gamma)) - ...
        xi*(x.*y) + H0); 
    Hx = @(x,y)(2*alpha*x - beta*(x./(sqrt(x.^2 + gamma))) - xi*y);
    Hy = @(x,y)(2*alpha*y - beta*(y./(sqrt(y.^2 + gamma))) - xi*x);
    
    e = (1/2)*(7/5)*((1 + Hx(x,y).^2).*(vx.^2) + ...
        (1 + Hy(x,y).^2).*(vy.^2) + ...
        2*(Hx(x,y).*Hy(x,y).*vx.*vy)) + g*H(x,y);
    
end
