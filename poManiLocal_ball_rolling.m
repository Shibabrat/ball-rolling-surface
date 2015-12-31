function [xW,x0W] = poManiLocal_ball_rolling(x0po,T,frac, ...
    stbl,dir,del,param,tmfd,n_mfd_traj) 

% [xW,x0W] = poManiLocal_ball_rolling(x0po,T,frac,stbl,dir,del,param,tmfd,n_mfd_traj)
%
% Compute points on (un)stable manifold for a periodic orbit (of period T)
% from eigenvalues of the monodromy PHI(T,0) matrix (denoted here as phi_T)
%
% output: 
% xW  = all computed points on the desired manifold
%
% x0W = initial points on stable or unstable manifold of a point
%	 	on the periodic orbit (used for globalizing)
%
% input:  
% x0po = initial (reference) point on the 3D periodic orbit in
% nondim. coords. at time t=0,T,2*T, etc.
%
% T    = period of periodic orbit in nondim. time 
% frac = {0 to 1}, fraction along orbit orbit beyond the reference
% point (point of which manifold is obtained) 
%
% stbl = +1 for unstable manifold
%	   = -1 for   stable manifold
%
% dir  = +1 for positive branch of manifold 
%	   = -1 for negative branch of manifold 
%
% del  = initial displacement along eigendirection (this varies with the
% size of the periodic orbit)
%
% dependencies:
%   [x,t,phi_tf,PHI] = stateTransitionMatrix_boatPR(x0,tf,R,OPTIONS,fixed_step) 
%   [x,t] = trajGet_boatPR(x0,tb,tf,R,OPTIONS)
%-----------------------------------------------------------------------
%
% Shane Ross (revised 2009-July-07), Shibabrat Naik (adapted and modified
% for roll-pitch ship model: 2014-May-20)

    if nargin < 9,
        n_mfd_traj = 10000 ;
        if nargin < 8,
            tmfd=3.7*T;
        end
    end

    N = 4; % dimension of phase space

    OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy
    % OPTIONS = odeset('RelTol',3e-10,'AbsTol',1e-12);  % lower accuracy

    % Get monodromy matrix, change function here 
    [xx,tt,phi_T,PHI] = stateTransitionMatrix_ball_rolling(x0po,T,OPTIONS,n_mfd_traj);

    n_mfd_traj_max = length(tt) ;

    if n_mfd_traj==10000,
        n_mfd_traj=n_mfd_traj_max;
    end


    % get eigenvalues of monodromy matrix and stable/unstable directions

    [sn,un,cn,y1Ws,y1Wu,y1Wc]   = eigGet(phi_T,1);


    disp(sprintf(':: Local Dynamics Summary'));
%     disp(sprintf('Parameters (R = %0.5g)',param));
    disp(        '================================');
    disp(sprintf('Number of   stable directions: %d',length(sn)))
    disp(sprintf('Number of unstable directions: %d',length(un)))
    disp(sprintf('Number of   center directions: %d',length(cn)))
    disp(        ' ');
    disp(sprintf('Number of computed pnts on po: %d',n_mfd_traj_max))
    disp(        ' ');
    disp(sprintf('Number of trajectories on mfd: %d',n_mfd_traj))
    disp(        ' ');


    % get purely real components of stable/unstable directions 

    numr=0;
    for k=1:length(sn), 
      if isreal(sn(k)), 
        numr=numr+1;
        snreal(numr)=sn(k);
        y1Wsreal(1:N,numr)=y1Ws(1:N,k);
      end
    end

    numr=0;
    for k=1:length(un), 
      if isreal(un(k)), 
        numr=numr+1;
        unreal(numr)=un(k);
        y1Wureal(1:N,numr)=y1Wu(1:N,k);
      end
    end

    % some 3D periodic orbits have 2-dim real stable and 
    % unstable manifolds (here we choose only one eigendirection)
    NN=1;
    WS = y1Wsreal(:,NN);  %   stable direction at initial point
    WU = y1Wureal(:,NN);  % unstable direction at initial point

    % now get the direction of desired manifold 
    mfrac=toTime(tt,frac*T);
    phi_frac=reshape(PHI(mfrac,1:N^2),N,N);

    % Which manifold, stable or unstable?
    if     stbl == -1,  % stable manifold Ws
      MAN = dir*phi_frac*WS;
      tb=-tmfd; tf=0;   % integrate backward
      colormfd='g';     % color green

    elseif stbl == +1,  % unstable manifold Wu
      MAN = dir*phi_frac*WU;
      tb=0; tf=tmfd;    % integrate forward
      colormfd='r';     % color red
    end


    % magnitude of displacement from point on periodic
    % as suggested by Barden[1994], (200km for Sun-Earth system)

    mag = norm(MAN(1:N/2)); 
    d = del/mag;	


    % point at t=frac*T (beyond reference pnt) on periodic orbit

    po_frac = xx(mfrac,:)'; 


    % point on desired manifold of p.o.

    x0W = po_frac + d*MAN; 


    % set accuracy for integration of individual trajectories on tube

    %OPTIONS = odeset('RelTol',3e-6,'AbsTol',1e-6); % lowest accuracy
    %OPTIONS = odeset('RelTol',3e-8,'AbsTol',1e-8); % med accuracy
    OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14); % high accuracy


    % integrate and plot the individual trajectories

    mfd_traj=0; 
    xW=[];
    xeU1 = [];
    numPts = 0;
    numPtsTraj = []; 
    % =======
    for mfrac = 1:round(n_mfd_traj_max/(n_mfd_traj-1)):n_mfd_traj_max,
        %  mfrac
        mfd_traj=mfd_traj + 1 ;
        phi_frac=reshape(PHI(mfrac,1:N^2),N,N);
        if      stbl == -1, % stable mfd
            MAN   = dir*phi_frac*WS ;  
        elseif  stbl == +1, % unstable mfd
            MAN   = dir*phi_frac*WU ;
        end
        mag = norm(MAN(1:N/2));       
        d = del/mag;

        po_frac = xx(mfrac,:)'; % point on p.o.
        x0W_all(mfrac,1:N) = po_frac(:)' + d*MAN(:)'; % displace onto manifold

        % Get trajectory on the tube manifolds
%         [x,t] = trajGet_boatPR(x0W_all(mfrac,1:N)',tb,tf,param,OPTIONS);
%         plot3(x(:,1), x(:,2), x(:,4), ['-',colormfd]);
% %         plot( x(end,2), x(end,4), 'bo'); 
%         if mfd_traj==1, 
%            hold on;
%            xlabel('$x$','interpreter','latex','fontsize',20);
%            ylabel('$y$','interpreter','latex','fontsize',20);
%         end
%         x0W(mfd_traj,1:N)=x0W_all(mfrac,1:N) ;
%         xW = [xW; x];
%         save mfd.dat xW x0W -ASCII -double
        
        
        % Get tube intersection with the Surface-Of-Section; takes the same
        % inputs as the above function, except also returns the event
        % solutions
        eventSwitch = 'on';
%         numTimesteps = 100;
%         [x,t,te,xe,ie] = get_traj_sos_boatPR(x0W_all(mfrac,1:N)',tb,tf,param,OPTIONS,eventSwitch,numTimesteps);
        [x,t,te,xe,ie] = get_traj_sos_ball_rolling(x0W_all(mfrac,1:N)',tb,tf,OPTIONS,eventSwitch);
%         plot( x(:,1), x(:,2),colormfd);
        plot3(x(:,1), x(:,2), x(:,4), ['-',colormfd]);
        
        if mfd_traj==1, 
           hold on;grid on;
           xlabel('$x$','interpreter','latex','fontsize',20);
           ylabel('$y$','interpreter','latex','fontsize',20);
           zlabel('$v_y$','interpreter','latex','fontsize',20);
        end
        x0W(mfd_traj,1:N)=x0W_all(mfrac,1:N) ;
        xW = [xW; x];
        numPtsTraj = [numPtsTraj; size(x,1)];
        numPts = numPts + 1;
        
        % Storing the event of intersection with plane U1
        if ~isempty(xe)
            xeU1 = [xeU1; xe(end,:)];
%             if (size(xe,1) >= 2)
%                 xeU1 = [xeU1; xe(3,:)]; 
%             end
        end

    end
    
%     tubeXU1 = permute(xeU1,[3 2 1]);
    save mfd.txt xW x0W -ASCII -double
    save xeU1.txt xeU1 -ASCII -double
    save mfd_params.txt numPtsTraj mfd_traj -ASCII -double
     
end
