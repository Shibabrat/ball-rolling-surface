function [x0po,T] = poBracketEnergy_ball_rolling(energyTarget,x0podata) 

%        [x0po,T] = poBracketEnergy_ball_rolling(energyTarget) ;
%
% Generate a family of periodic orbits (po) given a pair of seed initial
% conditions from a data file, while targeting a specific periodic orbit.
% This is performed using a scaling factor of the numerical continuation
% step of the phase space coordinates.
%
% Shane Ross (revised 2.19.04)
% Shibabrat Naik (modified for ball rolling on the surface: 26-Dec-2015)

    % delt = guessed change in period between successive orbits in family
    delt = -1.e-9 ;   % <==== may need to be changed
    energyTol = 1e-6;
    
    N = 4 ; % dimension of phase space

    %Load the first two members of the family from the data file
    if nargin < 2
        disp('Loading the periodic orbit family from data file')
        x0podata = importdata('x0po_T_energy_case4.txt');
    end
        
    x0po(1:2,1:N)   = x0podata(end-1:end,1:N);
    T(1:2,1)        = x0podata(end-1:end,N+1);
    energyPO(1:2,1) = get_energy_points_ball_rolling(x0po(1:2,:));
    iFam = 3;
    scaleFactor = 1;   %<=====scaling the change in initial guess
    finished = 1;
    
    while finished ~= 0, 
        
        FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
        disp(FAMNUM) ;
        
        %change in initial guess for next step
        dx  = x0po(iFam-1,1) - x0po(iFam-2,1) ;
        dy  = x0po(iFam-1,2) - x0po(iFam-2,2) ;

        %check p.o. energy and set the initial guess with scale factor
        if energyPO(iFam-1,1) < energyTarget,
            scaleFactor = scaleFactor
            x0po_g = [ (x0po(iFam-1,1) + scaleFactor*dx) ...
                (x0po(iFam-1,2) + scaleFactor*dy) 0 0] ;
            
            [x0po_iFam,tfpo_iFam] = poDifCor_ball_rolling(x0po_g) ; 
            energyPO(iFam,1) = get_energy_points_ball_rolling(x0po_iFam) ; 
            x0po(iFam,1:N) = x0po_iFam ;
            T(iFam,1)      = 2*tfpo_iFam ;
            iFam = iFam + 1;
            
            %to improve speed of stepping, increase scale factor for very
            %close p.o., this is when target p.o. hasn't been crossed,
            if abs(dx) < 1e-4 && abs(dy) < 1e-4,
                scaleFactor = 10;
            else
                scaleFactor = 1;
            end
            
        elseif energyPO(iFam-1,1) > energyTarget,
            scaleFactor = scaleFactor*1e-2
            x0po_g = [ (x0po(iFam-2,1) + scaleFactor*dx) ...
                (x0po(iFam-2,2) + scaleFactor*dy) 0 0] ;
            
            [x0po_iFam,tfpo_iFam] = poDifCor_ball_rolling(x0po_g) ; 
            energyPO(iFam-1,1) = get_energy_points_ball_rolling(x0po_iFam) 
            x0po(iFam-1,1:N) = x0po_iFam ;
            T(iFam-1,1)      = 2*tfpo_iFam ;
        
        end
        
        if abs(energyTarget - energyPO(iFam-1,1)) > ...
            energyTol*energyPO(iFam-1,1) + energyTol,
            finished = 1;
        else
            finished = 0;
        end
                    
    end

    POENERGYERR = sprintf('Error in the po energy from target %e', ...
        abs(energyTarget - energyPO(iFam-1,1))/energyPO(iFam-1,1)) ;
    disp(POENERGYERR) ;
        
    dum = [x0po T energyPO] ;
    save x0po_T_energy.txt -ascii -double dum

end
