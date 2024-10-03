function dUdt = rhs_lae(U, Pinv, D, E, a, sigma, g, t)
%
% dUdt = rhs_lae(U, Pinv, D, E, a, sigma, g, t)
% 
%  Implement right-hand-side computation for the linear advection equation
%  using an SBP-SAT scheme.
% 
%    INPUT: U - solution vector
%           Pinv - inverse of the integration matrix P
%           D - SBP differentiation matrix
%           E - convenience matrix for the SAT, e.g., E = diag(1, 0, ... ,0)
%           a - wave speed of the linear advection equation
%           sigma - penalty term of the SAT which must be <= -0.5
%           g(t) - (possibly) time dependent boundary data function
%           t - current time
% 
%    OUTPUT: dUdt - RHS of the SBP-SAT linear advection scheme

   % declare memory for output variable
   N = length(U);
   dUdt = zeros(N, 1);

   % compute the SAT term for a particular sigma
   SAT = a * sigma * Pinv * E * (U - [g(t); zeros(N-1,1)]);

   % assemble the right-hand-side
   % this approximates u_t = - a u_x + a sigma Pinv E (u_1 - g(t))
   dUdt = -a * D * U + SAT;

end