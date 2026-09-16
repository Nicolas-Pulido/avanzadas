clear;
clc;
close all;


Nx = 128;
Ny = 128;

Lx = 2*pi;
Ly = 2*pi;

hx = Lx/Nx;
hy = Ly/Ny;


x = linspace(-0.5*Lx + hx, 0.5*Lx, Nx);
y = linspace(-0.5*Ly + hy, 0.5*Ly, Ny);

[xx, yy] = ndgrid(x, y);



epsilon = 0.05;
Cahn = epsilon^2;


theta = atan2(yy, xx);
r = sqrt(xx.^2 + yy.^2);



r0 = 1.8 + 0.5*cos(4*theta) + 0.3*cos(8*theta);



u = tanh((r0 - r)/(sqrt(2)*epsilon));



p = 2*pi/Lx * [0:Nx/2-1, -Nx/2:-1];
q = 2*pi/Ly * [0:Ny/2-1, -Ny/2:-1];

p2 = p.^2;
q2 = q.^2;

[pp2, qq2] = ndgrid(p2, q2);



dt = 0.01;
T = 3;

Nt = round(T/dt);
ns = Nt/20;



figure(1);
clf;

contourf(x, y, real(u'), [0 0]);

axis image;
axis([x(1) x(Nx) y(1) y(Ny)]);

title('Condicion inicial modificada');
xlabel('x');
ylabel('y');

pause(0.01);



for iter = 1:Nt

    u = real(u);

    % Transformada de Fourier y termino no lineal
    s_hat = fft2(Cahn*u - dt*(u.^3 - 3*u));

    % Actualizacion en el dominio espectral
    v_hat = s_hat ./ ...
        (Cahn + dt*(2 + Cahn*(pp2 + qq2)));

    % Regreso al dominio espacial
    u = ifft2(v_hat);


    if mod(iter, ns) == 0

        contourf(x, y, real(u'), [0 0]);

        axis image;
        axis([x(1) x(Nx) y(1) y(Ny)]);

        title(['Tiempo = ', num2str(iter*dt)]);
        xlabel('x');
        ylabel('y');

        pause(0.01);

    end

end
