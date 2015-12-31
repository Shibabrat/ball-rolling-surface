function [x0po,tfpo] = po_auto_shooting_ball_rolling(x0poGuess,tfpoGuess)

%Defining the vector field
%constant parameters for the rolling surface    
alpha = 0.07;
beta = 1.017;
gamma = 15.103;
xi = 0.00656;
H0 = 12.065;
g = 918;

%function for the rolling surface 
% H = @(x,y)(alpha*(x^2 + y^2) - ...
%     beta*(sqrt(x^2 + gamma) + sqrt(y^2 + gamma)) - ...
%     xi*x*y + H0); 
Hx = @(x,y)(2*alpha*x - beta*(x/(sqrt(x^2 + gamma))) - xi*y);
Hy = @(x,y)(2*alpha*y - beta*(y/(sqrt(y^2 + gamma))) - xi*x);
Hxx = @(x,y)(2*alpha - (beta*gamma)/((x^2 + gamma)^(3/2)) );
Hyy = @(x,y)(2*alpha - (beta*gamma)/((y^2 + gamma)^(3/2)) );
Hxy = @(x,y)(-xi);

h = @(x,y,vx,vy)( (5*g)/7 + Hxx(x,y)*vx^2 + ...
    2*Hxy(x,y)*vx*vy + Hyy(x,y)*vy^2)/(1 + Hx(x,y)^2 + Hy(x,y)^2);


RelTol = 1e-12;     AbsTol = 1e-10;
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
% [x,t,phitf,PHI] = ...
%     stateTransitionMatrix_ball_rolling(x0po,2*tfpo,1,OPTIONS);
%   

x0po = x0poGuess;
tfpo = tfpoGuess;
iteration = 0;
iterMax = 50;
finished = 0;
Er = 1e-6; Ea = 1e-8;

while ~finished,

    if iteration > iterMax
        ERROR = 'Maximum iterations exceeded' ;
        disp(ERROR) ;
        break
    end

    [x,t,phitf,PHI] = ...
        stateTransitionMatrix_ball_rolling(x0po,tfpo,1,OPTIONS);

    rhstf = [x(end,3); ...
        x(end,4); ...
        -h(x(end,1),x(end,2),x(end,3),x(end,4))*Hx(x(end,1),x(end,2)); ...
        -h(x(end,1),x(end,2),x(end,3),x(end,4))*Hy(x(end,1),x(end,2))];

    rhst0 = [x(1,3); ...
        x(1,4); ...
        -h(x(1,1),x(1,2),x(1,3),x(1,4))*Hx(x(1,1),x(1,2)); ...
        -h(x(1,1),x(1,2),x(1,3),x(1,4))*Hy(x(1,1),x(1,2))];


    %linear system for updating the correction
    A = [phitf - eye(length(rhstf)) rhstf;
        rhst0' 0];

    b = [x(1,:)' - x(end,:)';
        0];

    deltaX = A\b;

    x0poNew = x0po + deltaX(1:4);
    tfpoNew = tfpo + deltaX(5);

    norm(deltaX(1:4))
    x0po = x0poNew;
    tfpo = tfpoNew;
    iteration = iteration + 1

    %Convergence test
    if (norm(deltaX(1:4)) < Er*norm(x0po) + Ea),
        finished = 1;
    end

end


end