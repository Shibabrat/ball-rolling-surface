%PLOT and TEST periodic orbit initial condition read from the ASCII data
%file
x0po_T = importdata('x0po_T_energy_case3.txt');
RelTol = 1e-12;     
AbsTol = 1e-14;
OPTIONS = odeset('RelTol',RelTol,'AbsTol',AbsTol); 
numTimesteps = 1e3;

%plot properties
axesFontName = 'factory';
% axesFontName = 'Times New Roman';
axFont = 18;
textFont = 18;
labelFont = 30;
lw = 2;    
set(0,'Defaulttextinterpreter','latex', ...
    'DefaultAxesFontName', axesFontName, ...
    'DefaultTextFontName', axesFontName, ...
    'DefaultAxesFontSize',axFont, ...
    'DefaultTextFontSize',textFont, ...
    'Defaultuicontrolfontweight','bold', ...
    'Defaulttextfontweight','bold', ...
    'Defaultaxesfontweight','bold');

%compute the periodic orbit 
po = struct();
for ii = 1:size(x0po_T,1)
    x0po = x0po_T(ii,1:4);
    tfpo = x0po_T(ii,5);
    e(ii) = get_energy_points_ball_rolling(x0po);
    
    [po(ii).tOut,po(ii).xOut] = ...
        ode113('ball_rolling2dof',linspace(0,tfpo,numTimesteps)',x0po,OPTIONS);
    
end
%%
figure(1)
for ii = 1:size(x0po_T)
    plot3(po(ii).xOut(:,1),po(ii).xOut(:,2),po(ii).xOut(:,3),'.-')
    grid on;hold on
    plot3(po(ii).xOut(1,1),po(ii).xOut(1,2),po(ii).xOut(1,3),'*b')
    plot3(po(ii).xOut(end,1),po(ii).xOut(end,2),po(ii).xOut(end,3),'ob')
    title(['Energy of the periodic orbit ',num2str(e(ii)),'(cm/s)$^2$']);
    xlabel('$x$');
    ylabel('$y$');
    zlabel('$v_x$');
    pause(1)
end
