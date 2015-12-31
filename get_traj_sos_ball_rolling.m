function [x,t,te,xe,ie] = get_traj_sos_ball_rolling(x0,tb,tf,OPTIONS,eventSwitch,numSteps)

    x = [];
    if nargin < 6
        TSPANtb = [-1.e-12 tb ] ;
        TSPANtf = [ 0      tf ] ;
    else
        TSPANtb = [-1.e-12:(tb+1e-12)/numSteps:tb] ;
        TSPANtf = [0:tf/numSteps:tf] ;
    end
    
    switch lower(eventSwitch)
        case 'off'
            MODEL = 'ball_rolling2dof_sos';

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
            % Event not tracked
            te = [];
            xe = [];
            ie = [];
        case 'on'
            OPTIONS = odeset('RelTol',1e-12,'AbsTol',1e-14,'Events','on'); 
            MODEL = 'ball_rolling2dof_sos';

            if tb==0,
              [t,x,te,xe,ie]     = ode113(MODEL,TSPANtf,x0,OPTIONS);
            elseif tf==0,
              [t,x,te,xe,ie]     = ode113(MODEL,TSPANtb,x0,OPTIONS);
            else
              [tt1,xx1,te1,xe1,ie1] = ode113(MODEL,TSPANtb,x0,OPTIONS);
              [tt2,xx2,te2,xe2,ie2] = ode113(MODEL,TSPANtf,x0,OPTIONS);

              x=[flipud(xx1);xx2];
              t=[flipud(tt1);tt2];
              te = [te1;te2];
              xe = [xe1;xe2];
              ie = [ie1;ie2];
            end
    end
    
end
