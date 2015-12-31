function plot_hills_region_rolling_ball(energyVal)

%PLOT_HILLS_REGION_ROLLING_BALL to plot the Hills region for a given energy
%for the ball rolling on a surface
%   plot_hills_region_rolling_ball(energyVal)
%
%   input:
%       energyVal: is the energy of the system and the Hills region plotted
%       corresponds to the realm of possible motion for this energy. 
 
    axesFontName = 'factory';
%     axesFontName = 'Times New Roman';
    set(0,'Defaulttextinterpreter','latex', ...
        'DefaultAxesFontName', axesFontName, ...
        'DefaultTextFontName', axesFontName, ...
        'DefaultAxesFontSize',24, ...
        'DefaultTextFontSize',24, ...
        'Defaultuicontrolfontweight','bold', ...
        'Defaulttextfontweight','bold', ...
        'Defaultaxesfontweight','bold');

    
    %constant parameters for the surface 
    alphaP = 0.07;
    betaP = 1.017;
    gammaP = 15.103;
    xiP = 0.00656;
    H0 = 12.065;    %cm
    g = 981;        %cm/s^2
    
    xMin = -15;
    xMax = 15;
    yMin = -15;
    yMax = 15;
    
    xPts = linspace(xMin,xMax,500)';
    yPts = linspace(yMin,yMax,500)';
    [xMesh, yMesh] = meshgrid(xPts, yPts);

    eMesh = g*(alphaP*(xMesh.^2 + yMesh.^2) - ...
        betaP*(sqrt(xMesh.^2 + gammaP) + sqrt(yMesh.^2 + gammaP)) - ...
        xiP*xMesh.*yMesh + H0);

    for i = 1:size(eMesh,1)
        for j = 1:size(eMesh,2)
            if (eMesh(i,j) > energyVal)
                colorMapHillsReg(i,j,:) = [0.5 0.5 0.5];
            else
                colorMapHillsReg(i,j,:) = [1 1 1];
            end
        end
    end

    mesh(xMesh, yMesh, eMesh/max(max(eMesh)), colorMapHillsReg)
    
    
    shading flat
    view(0,90)
    xlabel('$x$','Interpreter','Latex');
    ylabel('$y$','Interpreter','Latex','Rotation',0)
    axis equal tight
    box on
    
    set(gcf, 'PaperUnits', 'normalized')
    set(gcf, 'PaperPosition', [0 0 1 1])
    
    
end
