function Df = func_get_Jacobian(xEq)

    syms x y vx vy %h(x,y,vx,vy) H(x,y) Hx(x,y) Hxx(x,y) Jacobian(x,y,vx,vy)
  
    %constant parameters for the rolling surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
    H0 = 12.065;
    g = 918;
    xe = xEq(1);
    ye = xEq(2);
    vxe = xEq(3);
    vye = xEq(4);
    
    %function for the rolling surface 
    H(x,y) = alpha*(x^2 + y^2) - ...
        beta*(sqrt(x^2 + gamma) + sqrt(y^2 + gamma)) - ...
        xi*x*y + H0;
    Hx(x,y) = diff(H, x, 1);
    Hxx(x,y) = diff(Hx, x, 1);
    Hxy(x,y) = diff(Hx, y, 1);
    Hy(x,y) = diff(H, y, 1);
    Hyy(x,y) = diff(Hy, y, 1);
       
    h(x,y,vx,vy) = ( (5*g)/7 + Hxx(x,y)*vx^2 + ...
        2*Hxy(x,y)*vx*vy + Hyy(x,y)*vy^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2);
    hx(x,y,vx,vy) = diff(h, x, 1);
    hy(x,y,vx,vy) = diff(h, y, 1);
    hvx(x,y,vx,vy) = diff(h, vx, 1);
    hvy(x,y,vx,vy) = diff(h, vy, 1);

    %Computing the entries of the first derivative matrix at (xe,ye,0,0)
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

    %Use vector differentiation 
%     f(x,y,vx,vy) = [vx; vy; -h(x,y,vx,vy)*Hx(x,y); -h(x,y,vx,vy)*Hy(x,y)];
%     DfMat(x,y,vx,vy) = jacobian(f,[x y vx vy]);
%     Df = double(DfMat(xe,ye,vxe,vye));
    
end