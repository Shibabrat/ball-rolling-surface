%SCRIPT to generate initial set on the Poincaré S-O-S for safe set
%sculpting

%constant parameters for the rolling surface 
alpha = 0.07;
beta = 1.017;
gamma = 15.103;
xi = 0.00656;
H0 = 12.065;
g = 981;
E = 3304.24 + 30;

%function for the rolling surface 
H = @(x,y)(alpha*(x.^2 + y.^2) - ...
    beta*(sqrt(x.^2 + gamma) + sqrt(y.^2 + gamma)) - ...
    xi*(x.*y) + H0); 
Hx = @(x,y)(2*alpha*x - beta*(x./(sqrt(x.^2 + gamma))) - xi*y);
Hy = @(x,y)(2*alpha*y - beta*(y./(sqrt(y.^2 + gamma))) - xi*x);

%Define the points on the Poincaé S-O-S
yDomain = linspace(0,10,200);
vxDomain = linspace(-40,40,200);
xDomain = yDomain;
[xMesh,yMesh,vxMesh] = ndgrid(xDomain,yDomain,vxDomain);

%Solve for vy-coordinates using the quadratic
aQuad = (1/2)*(7/5)*(1 + Hx(xMesh,yMesh).^2);
bQuad = (1/2)*(7/5)*(2*Hx(xMesh,yMesh).^2.*vxMesh);
cQuad = g*H(xMesh,yMesh) - E + ...
    (1/2)*(7/5)*(1 + Hx(xMesh,yMesh).^2).*(vxMesh.^2);
vyMeshP = (-bQuad + sqrt(bQuad.^2 - 4*(aQuad.*cQuad)))./(2*aQuad);
vyMeshM = (-bQuad - sqrt(bQuad.^2 - 4*(aQuad.*cQuad)))./(2*aQuad);

%Keep the vy-coordinates that are real and satisfy constraint to lie on the
%S-O-S
ll = 1;
for ii = 1:size(vyMeshP)
    for jj = 1:size(vyMeshP)
        for kk = 1:size(vyMeshP)
            if (imag(vyMeshP(ii,jj,kk)) > 1e-12) && ...
                    (real(vyMeshP(ii,jj,kk)) > vxMesh(ii,jj,kk))
                gridPts(ll,1:2) = [xMesh(ii,jj,kk), vxMesh(ii,jj,kk)];
                ll = ll + 1;
            else
                ll = ll;
            end
        end
    end
end




