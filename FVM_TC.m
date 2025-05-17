%% -----PENALIZATION METHOD FOR 1D DIFFUSION EQUATION----- %%

clear all;
close all;

% Grid elements
Nx = 512;
L = 4.0;
dx = L/Nx;
xf = linspace(-2,2,Nx+1);
xc = 0.5*(xf(1:end-1)+xf(2:end));

% Model constants
k = 1;
eta = 1e-2;

% Time parameter
dt = 1e-4;
tmax = 1;
Nt = round(tmax/dt);
si = 100;

% Mask
xi = ones(1,Nx);
for i=1:Nx
    if (xc(i) >= -1) && (xc(i) <= 1)
        xi(i) = 0;
    end
end

[~,idn1] = min(abs(xf - (-1)));
[~,idn2] = min(abs(xf - (1)));

tht = cos(4*pi*xc) + cos(pi*xc);
tht0 = tht;
Df = zeros(1,Nx+1);
Fx = zeros(1,Nx+1);

t = 0;
fprintf('Starting time integration...\n');

for n=1:Nt

    for i=2:Nx
        Dl = eta*xi(i-1) + k*(1-xi(i-1));
        Dr = eta*xi(i) + k*(1-xi(i));

        Df(i) = 0.5*(Dl+Dr);
    end
    
    Dl = eta*xi(Nx) + k*(1-xi(Nx));
    Dr = eta*xi(1) + k*(1-xi(1));

    Df(1) = 0.5*(Dl + Dr);
    Df(Nx+1) = Df(1);

    Df(idn1) = 0;
    Df(idn2) = 0;

    diag = zeros(Nx,3);

    for i=1:Nx
        diag(i,2) = -(Df(i)+Df(i+1))/dx^2;
    end

    for i=2:Nx
        diag(i,1) = Df(i)/dx^2;
    end
    diag(1,1) = Df(1)/dx^2;

    for i=1:Nx-1
        diag(i,3) = Df(i+1)/dx^2;
    end
    diag(Nx,3) = Df(Nx+1)/dx^2;

    e = ones(Nx,1);
    L = spdiags([diag(:,1), diag(:,2), diag(:,3)], [-1, 0, 1], Nx, Nx);

    L(1, Nx) = Df(1)/dx^2;
    L(Nx, 1) = Df(Nx+1)/dx^2;

    exp_l = expm(full(L)*dt);
    tht = exp_l * tht';
    tht = tht';
    thte = exp(-16*pi^2*k*t)*cos(4*pi*xc) + exp(-pi^2*k*t)*cos(pi*xc);

    t = t + dt;    
    % Plot and save results at specified intervals
    if mod(n, si) == 0 || n == Nt
        fprintf('Time step %d of %d (t = %.3f)\n', n, Nt, t);
        
            hfig = figure('Visible','off');
            hold on; grid on;
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
            ylabel('$\theta_{\eta}$', 'Interpreter', 'latex', 'FontSize', 15);
            
            plot(xc, tht, '--', 'LineWidth', 1.5);
            plot(xc, thte, 'Color', [0.4,0.0,0.2], 'LineWidth', 1.5);
            ylim([-2,2]);
            title(sprintf('t = %.3f', t));
            filename = sprintf('frame_%05d.png', n);
            saveas(hfig, filename);
            close(hfig); 
    end
end

%% Solution comparison

figure;
plot(xc, tht0, color=[0.0,0.2,0.6], LineWidth=1.5);
hold on;
plot(xc, tht, color=[0.6,0.0,0.6], LineWidth= 1.5);
xline(-1,'--');
xline(1,'--');
ylim([-2,2]);
xlabel('$x$',Interpreter='latex',FontSize=15);
ylabel('$\theta_{\eta}$',Interpreter='latex',FontSize=15);
legend('Initial condition', 'Final solution','Neumann boundary','Neumann boundary');
grid on;

% Mask function

figure;
plot(xc, xi, 'k-', 'LineWidth', 2);
grid on;
xlabel('$x$',Interpreter='latex',FontSize=15);
ylabel('$\chi$',Interpreter='latex',FontSize=15);
title('Mask function (1 = solid, 0 = fluid)');
ylim([-0.1, 1.1]);

fprintf('Simulation completed.\n');