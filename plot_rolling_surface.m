function plot_rolling_surface

%PLOT_ROLLING_SURFACE plots a surface for the ball rolling problem 

    axesFontName = 'factory';
%     axesFontName = 'Times New Roman';
    set(0,'Defaulttextinterpreter','latex', ...
        'DefaultAxesFontName', axesFontName, ...
        'DefaultTextFontName', axesFontName, ...
        'DefaultAxesFontSize',20, ...
        'DefaultTextFontSize',30, ...
        'Defaultuicontrolfontweight','bold', ...
        'Defaulttextfontweight','bold', ...
        'Defaultaxesfontweight','bold');

    
    %constant parameters for the surface 
    alpha = 0.07;
    beta = 1.017;
    gamma = 15.103;
    xi = 0.00656;
    H0 = 12.065;

    xMin = -15;
    xMax = 15;
    yMin = -15;
    yMax = 15;
    
    xPts = linspace(xMin,xMax,100)';
    yPts = linspace(yMin,yMax,100)';
    
    [xMesh,yMesh] = meshgrid(xPts,yPts);
    
    
    rollingSurf = alpha*(xMesh.^2 + yMesh.^2) - ...
        beta*(sqrt(xMesh.^2 + gamma) + sqrt(yMesh.^2 + gamma)) - ...
        xi*xMesh.*yMesh + H0;
    
    surf(yMesh,xMesh,rollingSurf)
    xlabel('$x (cm)$')
    ylabel('$y (cm)$')
    zlabel('$H(x,y) (cm)$')
    grid on 
    shading interp
    axis tight equal
    colorbar
    
    set(gcf, 'PaperUnits', 'normalized')
    set(gcf, 'PaperPosition', [0 0 1 1])

end