function [x0po,T] = poFamGet_ball_rolling(eqNum,Ax1,Ax2,nFam) 

%        [x0po,T] = poFamGet_ball_rolling(eqNum,Ax1,Ax2,nFam) ;
%
% Generate a family of periodic orbits (po) given a pair of 
% seed initial conditions and periods
%
% Shane Ross (revised 2.19.04)
% Shibabrat Naik (modified for ball rolling on the surface: 18-Dec-2015)

    % delt = guessed change in period between successive orbits in family

    delt = -1.e-9 ;   % <==== may need to be changed

    N = 4 ; % dimension of phase space

    x0po = zeros(nFam,N) ;
    T    = zeros(nFam,1) ;

    [x0poGuess1,TGuess1] = poGuessLinear_ball_rolling(eqNum,Ax1) 
    [x0poGuess2,TGuess2] = poGuessLinear_ball_rolling(eqNum,Ax2)

% Get the first two periodic orbit initial conditions
    iFam = 1 ;
    FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
    disp(FAMNUM) ;
  [x0po1,tfpo1] = poDifCor_ball_rolling(x0poGuess1) ; 
%     [x0po1,tfpo1] = po_auto_shooting_ball_rolling(x0poGuess1, TGuess1);

    iFam = 2 ;
    FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
    disp(FAMNUM) ;
    [x0po2,tfpo2] = poDifCor_ball_rolling(x0poGuess2) ;               
%     [x0po2,tfpo2] = po_auto_shooting_ball_rolling(x0poGuess2, TGuess2);

    x0po(1:2,1:N) = [x0po1(:)'  ; x0po2(:)' ];
    T   (1:2)     = [2*tfpo1      ; 2*tfpo2     ]; 

    %Generate the other members of the family using numerical continuation
    for iFam = 3:nFam,

        FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
        disp(FAMNUM) ;

        dx  = x0po(iFam-1,1) - x0po(iFam-2,1) ;
        dy  = x0po(iFam-1,2) - x0po(iFam-2,2) ;
%         dOrbit = x0po(iFam-1,:) - x0po(iFam-2,:);
        dt  = T(iFam-1) - T(iFam-2) ;

%         x0po_g = x0po(iFam-1,:)' + dOrbit';
        t1po_g =   (T(iFam-1) + dt)/2 + delt ;
        x0po_g = [ (x0po(iFam-1,1) + dx) (x0po(iFam-1,2) + dy) 0 0] ;
%         t1po_g = T(iFam-2) + abs(dt);

      % differential correction takes place in the following function
        [x0po_iFam,tfpo_iFam] = poDifCor_ball_rolling(x0po_g) ; 
%         [x0po_iFam,tfpo_iFam] = ...
%             po_auto_shooting_ball_rolling(x0po_g', t1po_g);

        x0po(iFam,1:N) = x0po_iFam ;
%         T   (iFam)     = 2*t1_iFam ;	
        T(iFam)        = 2*tfpo_iFam;

      if mod(iFam,10) == 0,
        dum = [x0po T] ;
    %     save x0po_T.dat -ascii -double dum
      end

    end

dum = [x0po T] ;
save x0po_T.txt -ascii -double dum

end
