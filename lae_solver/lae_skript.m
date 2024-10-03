
% driver script for approximating the linear advection equation
% solution via an SBP-SAT scheme on a single block
%
%    u_t + a u_x = 0, a is a positive real constant on [x_l,x_r]
%    u(x,0) = u_0(x) % initial condition
%    u(a,t) = g(t)   % boundary data

% set the domain
x_r = 3.0;
x_l = 0.0;

% resolution and get x points
N = 100;

% create uniform gridpoints
dx = (x_r - x_l) / (N - 1);
x = transpose(x_l:dx:x_r);

% create the SBP operator pair
% [P, D] = sbp42(N, dx);
[P, D] = sbp63(N, dx);
% [P, D] = sbp84(N, dx);

% store P inverse for convenience
Pinv = zeros(N,N);
for j = 1:N
   Pinv(j,j) = 1 / P(j,j);
end

% boundary matrix to build the SAT on the left with E=diag(1,0,...,0,0)
% OBS! if you switch a to be negative then the boundary data is set at the
%      right and this matrix would be E=diag(0,0,...,0,1)
E = zeros(N,N);
E(1,1) = 1;


% Set the SAT penalty parameter. For LAE must be <= -1/2
sigma = -1.0;

a = 1.1; % wavespeed

% % setup the manufatured solution and boundary term
% u_ex = @(x,t) 2 + sin(2 * pi * (x - a * t));
% g = @(t) 2 + sin(2 * pi * (x_l - a * t));

% setup the manufatured solution and boundary term
% set the domain to be [0,3]
u_ex = @(x,t) exp(-20*((x - a * t) - 1.5).^2);
g = @(t) 0.0;


% % Try something fancier that emphasizes what the boundary condition does
% u_ex = @(x,t) zeros(size(x));
% g = @(t) boundary_condition_sine_sector(t);

% set the initial condition
U = u_ex(x, 0.0);


% set the time step size
CFL = 1.0; % for SBP 42 or 63 or 84
dt = CFL * dx / abs(a);
t = 0.0;
t_final = 5.0; % 1.0 for error analysis
k = 0; % # of time steps (only used for plotting data)


% Do the time loop
while t < t_final
   % Avoid stepping over the final time because we use a while loop
   if t + dt > t_final
     dt = t_final - t;
   end
   t = t + dt;
   k = k+1;
   U = step_by_rk3(t, dt, U, Pinv, D, E, a, sigma, g);

   % plot on the fly to show a "movie"
   % OBS! somewhat ugly way to do this
   if mod(k, 5) == 0
      % plot(x,u_ex(x,t),'-k','LineWidth',1.5)
      plot(x, U, '-k', 'LineWidth', 1.5)
      hold on
      plot(x_l, g(t), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
      hold off
      xlim([x(1) x(end)])
      % ylim([1.0 3.0])
      ylim([-1.0 1.2])
      xlabel('$x$', 'interpreter', 'latex')
      % legend('$u^{\textrm{exact}}$', '$u$', 'interpreter', 'latex')
      legend('$U$', '$g(t)$', 'interpreter', 'latex')
      title(['$t = $',num2str(t)], 'interpreter', 'latex')
      set(gca, 'fontsize', 24)

      pause(0.005)
   end
end

%% Commands to check the operator spectrum

% draw the stability region for RK3 (or forward Euler)
[x, y] = meshgrid(-3:0.1:1, -3:0.1:3);
z = x + 1i * y;

% % forward Euler stability region
% FE = 1 + z; 
% FEModulus = abs(FE);
% contourf(x, y, -FEModulus, [-1 -1])

% Runge-Kutta 3 stability region
RK3 = 1 + z + z.^2/2 + z.^3/6;
RK3Modulus = abs(RK3);
contourf(x, y, -RK3Modulus, [-1 -1])

% setup nice axes and labels
hold on

plot([x(1) x(end)], [0 0], '-k')
plot([0 0], [y(1) y(end)], '-k')
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
set(gca, 'FontSize',16)
% axis equal

% Compute operator spectrum that collects everything that "hits" U
  operator = -a * D + a * sigma * Pinv * E;
  lamb = eig(operator);
% raw spectra
  % plot(real(lamb), imag(lamb), 'ro', 'MarkerFaceColor', 'r')
% spectra scaled by dt
  CFL = 1.0;
  dt = CFL * dx / a;
  plot(real(lamb)*dt, imag(lamb)*dt, 'ro', 'MarkerFaceColor', 'r')

hold off
