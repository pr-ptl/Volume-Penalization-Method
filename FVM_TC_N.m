%% -----PENALIZATION METHOD FOR 1D DIFFUSION EQUATION: ERROR VS N ----- %%

clear all; close all;

% Grid elements
N = round(2.^[3:0.5:15]);
err = zeros(size(N));

for ns=1:length(N)
    Nx = N(ns);
    L = 4.0;
    dx = L/Nx;
    xf = linspace(-2,2,Nx+1);
    xc = 0.5*(xf(1:end-1)+xf(2:end));

    % Model constants
    k = 1;

    % Time parameter
    dt = 1e-4;
    tmax = 1;
    Nt = round(tmax/dt);

    % Mask
    xi = ones(1,Nx);
    for i=1:Nx
        if (xc(i) >= -1) && (xc(i) <= 1)
            xi(i) = 0;
        end
    end

    [~,idn1] = min(abs(xf - (-1)));
    [~,idn2] = min(abs(xf - (1)));
    % cell points in fluid
    idxc = find(xc>=-1 & xc<=1);

    % Initial condition
    tht0 = cos(4*pi*xc) + cos(pi*xc);

    % Penalty parameters to test
    etav = 1e-8;

    fprintf('Running simulations for different grid size...\n');

    eta = etav;
    tht = tht0;
    t = 0;
    
    for n=1:Nt
        Df = zeros(1,Nx+1);
        for i=2:Nx
            Dl = eta*xi(i-1) + k*(1 - xi(i-1));
            Dr = eta*xi(i) + k*(1 - xi(i));
            Df(i) = 0.5*(Dl + Dr);
        end
        Dl = eta*xi(Nx) + k*(1 - xi(Nx));
        Dr = eta*xi(1) + k*(1 - xi(1));
        Df(1) = 0.5*(Dl + Dr);
        Df(Nx+1) = Df(1);

        Df(idn1) = 0;
        Df(idn2) = 0;

        diag = zeros(Nx,3);
        for i=1:Nx
            diag(i,2) = -(Df(i) + Df(i+1))/dx^2;
        end
        for i=2:Nx
            diag(i,1) = Df(i)/dx^2;
        end
        diag(1,1) = Df(1)/dx^2;

        for i=1:Nx-1
            diag(i,3) = Df(i+1)/dx^2;
        end
        diag(Nx,3) = Df(Nx+1)/dx^2;

        L = spdiags([diag(:,1), diag(:,2), diag(:,3)], [-1, 0, 1], Nx, Nx);
        L(1, Nx) = Df(1)/dx^2;
        L(Nx, 1) = Df(Nx+1)/dx^2;

        exp_l = expm(full(L)*dt);
        tht = exp_l * tht';
        tht = tht';
        t = t + dt;
    end
    
    thte = exp(-16*pi^2*k*t)*cos(4*pi*xc) + exp(-pi^2*k*t)*cos(pi*xc);
    err(ns) = sqrt(sum((tht(idxc) - thte(idxc)).^2) * dx);
    
    fprintf('N = %.1e, L2 error = %.3e\n', Nx, err(ns));
end

%% Plot error vs N

figure;
loglog(N, err, 'o-', 'LineWidth', 2);
xlabel('$N$',Interpreter='latex',FontSize=15);
ylabel('$||\epsilon||_2$',Interpreter='latex',FontSize=15);
title('Error vs Number of points $N$',Interpreter='latex',FontWeight='bold');
grid on;
exportgraphics(gcf,'l2normgrid.jpg','Resolution',500);