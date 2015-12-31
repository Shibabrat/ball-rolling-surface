function PHIdot = varEqns_ball_rolling(t,PHI)

%        PHIdot = varEqns_ball_rolling(t,PHI) ;
%
% This here is a preliminary state transition, PHI(t,t0),
% matrix equation attempt for a ball rolling on the surface, based on...
%
%        d PHI(t, t0)
%        ------------ =  Df(t) * PHI(t, t0)
%             dt
%
%-----------------------------------------------------------
% BALL ROLLING ON THE SURFACE CONVENTION:
%
%
%    2nd quadrant   L12(EQNUM=3)    1st quadrant
%
%
%     L23(EQNUM=4)                  L41(EQNUM=2)
%
%
%    3rd quadrant   L34(EQNUM=1)    4th quadrant
%-----------------------------------------------------------
% Shane Ross (revised 2.19.04)
% Shibabrat Naik (added pitch-roll ship model: 2014-May-20)

    R = 1.6;
    
% R: Ratio of pitch and roll frequencies
    
    x(1:4) = PHI(17:20);
    phi  = reshape(PHI(1:16),4,4);

    % The following is the Jacobian matrix
%     Df = func_get_Jacobian(x(1:4));

    %constant parameters for the rolling surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
    H0 = 12.065;
    g = 918;

    %function for the rolling surface 
    xe = x(1);
    ye = x(2);
    vxe = x(3);
    vye = x(4);
    H = @(x,y)(alpha*(x^2 + y^2) - ...
        beta*(sqrt(x^2 + gamma) + sqrt(y^2 + gamma)) - ...
        xi*x*y + H0); 
    Hx = @(x,y)(2*alpha*x - beta*(x/(sqrt(x^2 + gamma))) - xi*y);
    Hy = @(x,y)(2*alpha*y - beta*(y/(sqrt(y^2 + gamma))) - xi*x);
    Hxx = @(x,y)(2*alpha - (beta*gamma)/((x^2 + gamma)^(3/2)) );
    Hyy = @(x,y)(2*alpha - (beta*gamma)/((y^2 + gamma)^(3/2)) );
    Hxy = @(x,y)(-xi);
    Hyx = @(x,y)(-xi); 
    Hxxx = @(x,y)( (3*beta*gamma*x)/((x^2 + gamma)^(5/2)) );
    Hyyy = @(x,y)( (3*beta*gamma*y)/((y^2 + gamma)^(5/2)) );
    
    h = @(x,y,vx,vy)( (5*g)/7 + Hxx(x,y)*vx^2 + ...
        2*Hxy(x,y)*vx*vy + Hyy(x,y)*vy^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2);
    hx = @(x,y,vx,vy)( (Hxxx(x,y)*vx^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2) - ...
        h(x,y,vx,vy)*(2*Hx(x,y)*Hxx(x,y) + 2*Hy(x,y)*Hyx(x,y)) ); 
    hy = @(x,y,vx,vy)( (Hyyy(x,y)*vy^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2) - ...
        h(x,y,vx,vy)*(2*Hy(x,y)*Hyy(x,y) + 2*Hx(x,y)*Hxy(x,y)) ); 
    hvx = @(x,y,vx,vy)( (2*Hxx(x,y)*vx + 2*Hxy(x,y)*vy)/(1 + Hx(x,y)^2 + Hy(x,y)^2) );
    hvy = @(x,y,vx,vy)( (2*Hyy(x,y)*vy + 2*Hxy(x,y)*vx)/(1 + Hx(x,y)^2 + Hy(x,y)^2) );
   
    Df(1,:) = [0, 0, 1, 0];
    Df(2,:) = [0, 0, 0, 1];
    Df(3,1) = double(-hx(xe,ye,vxe,vye)*Hx(xe,ye) - h(xe,ye,vxe,vye)*Hxx(xe,ye));
    Df(3,2) = double(-hy(xe,ye,vxe,vye)*Hx(xe,ye) -h(xe,ye,vxe,vye)*Hxy(xe,ye));
    Df(3,3) = double(-hvx(xe,ye,vxe,vye)*Hx(xe,ye));
    Df(3,4) = double(-hvy(xe,ye,vxe,vye)*Hx(xe,ye));
    Df(4,1) = double(-hx(xe,ye,vxe,vye)*Hy(xe,ye) -h(xe,ye,vxe,vye)*Hxy(xe,ye));
    Df(4,2) = double(-hy(xe,ye,vxe,vye)*Hy(xe,ye) -h(xe,ye,vxe,vye)*Hyy(xe,ye));
    Df(4,3) = double(-hvx(xe,ye,vxe,vye)*Hy(xe,ye));
    Df(4,4) = double(-hvy(xe,ye,vxe,vye)*Hy(xe,ye));

%     Df    =[  0     0    1    0 ;
%               0     0    0    1 ;
%             (-1+2*x(2))  2*x(1)   0    0 ;
%             R^2*x(1)  -R^2  0    0 ];

    phidot = Df * phi; % variational equation

    PHIdot        = zeros(20,1);
    PHIdot(1:16)  = reshape(phidot,16,1); 
    PHIdot(17)    = x(3);
    PHIdot(18)    = x(4);
    PHIdot(19)    = -h(x(1),x(2),x(3),x(4))*Hx(x(1),x(2)) ;
    PHIdot(20)    = -h(x(1),x(2),x(3),x(4))*Hy(x(1),x(2));

end