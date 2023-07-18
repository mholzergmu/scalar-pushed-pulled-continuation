# scalar-pushed-pulled-continuation
Code to continue pulled fronts, pushed fronts and pushed to pulled transition points to accompany Avery, Holzer, Scheel, Pushed-to-pulled front transitions: continuation, speed scalings, and hidden monotonicty 2023. 

Performs pulled front continuation, pushed front continuation and locates pushed-to-pulled transition value for the scalar PDE 
%    u_t=u_{xx}+f(u,u_x,mu)

%%% Computes solutions to the BVP u_{xx}+cu_x+f(u,u_x,mu)=0 using the
%    decomposition u=  um*Chi_m+up*Chi_p+ Chi_p.*(a*x+b).*exp(nu*x)+ w

%   um -- stable homogeneous state in the wake
%   up -- unstable homogeneous state ahead of front 
%   (a*x+b).*exp(nu*x) -- far-field exponential ansatz
%   w -- localized core function 
%   Chi_m and Chi_p -- internally defined cut-off functions 

%%%%  Code requires five files
   %    DiffLap.m creates finite difference Laplacian and Derivative
   %    matricies

   %    scalar_cont.m is the main continuation code where all user
   %    specified changes take place

   %    scalar_sec.m specifies the nonlinear function F to be solved

   %    scalar_dF_sec.m computes the Jacobian 

   %    scalar_newton_sec.m facilitates Newton based continuation 

%%%% User must specify the following

    % the reaction function f(u,u_x,mu) and its derivatives f_u(u,u_x,mu)
    % and f_{u_x}(u,u_x,mu)

    % the dispersion relation for the unstable state d(lambda,nu) and its
    % nu derivative d_{nu}(\lambda,nu) in a co-moving frame with speed c

    % The default is to define L as the finite difference diffusion
    % operator with diffusion coefficient 1.  User has option to scale the
    % diffusion or include fourth order derivatives 

    % initialization of variables w, a, b, nu, c, up, um and mu

    % continuation parameters 

    % par - 1(pulled front continuation) , 2(pushed-to-pulled transition 
    %       location) or 4 (pushed front continuation) 

    % parseq - to perform a combination of continuation routines 

    % direction - whether to initially  continue in the positive (1) 
    %             or negative (-1) direction 


%%% Outputs 


%   mu1,mu2,mu4 -- mu parameter found during continuation for par=1,2 or 4
%   a1,a2, a4 -- a parameter found during continuation for par=1,2 or 4
%   b1,b2, b4 -- b parameter found during continuation for par=1,2 or 4
%   c1,c2, c4 -- wavespeed calculated during continuation for par=1,2 or 4
