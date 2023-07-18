clear all
close all

%%%%% Code to accompany Avery, Holzer and Scheel 2022 Pushed-to-pulled 
% front transitions: continuation, speed scalings, and hidden monotonicity

%%% Performs pulled front continuation, pushed front continuation and
%    locates pushed-to-pulled transition value for the scalar PDE 
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



format long;

%%% User inputs 


% Construct our space

M=20;
dx=.1;

N=floor(2*M/dx);
x=linspace(-M,M,N)';
dx=x(2)-x(1);

% define nonlinearity and dispersion relation 

%% Example: Nagumo 

f=@(u,u_x,mu) u.*(mu+u).*(1-u);
duf=@(u,u_x,mu) mu+2*u.*(1-mu)-3*u.^2; % derivative with respect to u
duxf=@(u,u_x,mu) 0*u; % derivative with respect to u_x

d= @(nu,c,mu) nu.^2+c*nu+mu; % dispersion relation 
dn=@(nu,c,mu) 2*nu+c;  % nu derivative of dispersion relation in co-moving frame 



%% Example: Burger-KPP

% f=@(u,u_x,mu) u.*(1-u)-mu*u.*u_x;
% duf=@(u,u_x,mu) 1-2*u-mu*u_x;
% duxf=@(u,u_x,mu) -mu*u;
% 
% d= @(nu,c,mu) nu.^2+c*nu+1;
% dn=@(nu,c,mu) 2*nu+c;

%% Example: EFKPP with Nagumo 

% f=@(u,u_x,mu) u+mu*u.^2-10*u.^3;
% duf=@(u,u_x,mu) 1+2*mu*u-30*u.^2; % derivative
% duxf=@(u,u_x,mu) 0*u;
% 
% zeta=1/75;
% d= @(nu,c,mu) -(zeta)*nu.^4+nu.^2+c*nu+1;
% dn=@(nu,c,mu) -(zeta)*4*nu.^3+2*nu+c;

% initialization 


%%% Nagumo

uinit=zeros(N+7,1); % w 
uinit(N+1)=1; %a
uinit(N+2)=.1; %b
uinit(N+3)=0; %up
uinit(N+4)=1; %um
uinit(N+5)=2; %c
uinit(N+6)=-1; %nu
uinit(N+7)=1; %mu

%% Bugers-KPP

% uinit=zeros(N+7,1);
% uinit(N+1)=1; %a
% uinit(N+2)=.1; %b
% uinit(N+3)=0; %up
% uinit(N+4)=1; %um
% uinit(N+5)=1; %c
% uinit(N+6)=-1; %nu
% uinit(N+7)=0; %mu

%% Extended Nagumo/Fisher-KPP

% uinit=zeros(N+7,1);
% uinit(N+1)=.16; %a
% uinit(N+2)=-0.01; %b
% uinit(N+3)=0; %up
% uinit(N+4)=sqrt(10); %um
% uinit(N+5)=1.9; %c
% uinit(N+6)=-1.12; %nu
% uinit(N+7)=0.001; %mu


% Create Laplacian and Derivative Matricies

fourth_order=true;
neumann=true; 

[D1,D2]=DiffLap(N,fourth_order,neumann,dx);

%L=D2-(zeta)*D2*D2; 
L=D2;

% Continuation 

% par==1 is pulled front continuation in parameter mu 
% par==2 is one step of pulled front with a=0 to locate transition point
% par==4 is pushed front continuation 

%parseq=[1 4]; % continue pulled front first then switch to follow pushed fronts
 parseq=[1 2]; % continue pulled fronts in mu to find transition point
%parseq=4; continue pushed fronts directly 

% numerical parameters
tol = 1e-10 ; 
max_it = 1000;
verbose=true;
solution=false;
direction=-1; %1 for increasing mu and -1 for decreasing mu 
tang_init=zeros(N+7,1);
tang_init(N+7)=direction;




% define cutoff functions
or=2; % for fourth order scheme, 1 enough for first order, should not hurt if too large
int0=[zeros(or,1);ones(N-2*or,1);zeros(or,1)]; % only 2 
eta=5;
%Chi+
Chi_p=1./(1+exp(-eta*x));
D1_Chi_p=int0.*(D1*Chi_p);%(eta*exp(eta*(-L/2-x)))./((1+exp(eta*(-L/2-x))).^2);
L_Chi_p=int0.*(L*Chi_p);% in case of D2 (eta^2*exp(eta*(-L/2-x))).*(exp(eta*(-L/2-x))-1)./((1+exp(eta*(-L/2-x))).^3);

Chi_m=1./(1+exp(eta*x));
D1_Chi_m=int0.*(D1*Chi_m);%(eta*exp(eta*(-L/2-x)))./((1+exp(eta*(-L/2-x))).^2);
L_Chi_m=int0.*(L*Chi_m);% in case of D2 (eta^2*exp(eta*(-L/2-x))).*(exp(eta*(-L/2-x))-1)./((1+exp(eta*(-L/2-x))).^3);

% define commutators
commLChip=spdiags(int0,0,N,N)*(L*spdiags(Chi_p,0,N,N)-spdiags(Chi_p,0,N,N)*L);


sec=zeros(N+7,1);
sec(N+7)=1;







% Initialize secant continuation


par=parseq(1);

U_o=uinit;



[U_o,status,steps]=...
scalar_newton_sec(U_o,D1,L,N,sec,tol,max_it,x,verbose,solution,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);


%return

status



figure(3)
U_n=U_o;


  up=U_n(N+3);
   um=U_n(N+4);
   c=U_n(N+5);
   nu=U_n(N+6);


init_guess_u=um*Chi_m+up*Chi_p+ Chi_p.*(U_n(N+1)*x+U_n(N+2)).*exp(nu*x)+ U_n(1:N);



  



plot(x,init_guess_u,'b');






[U_n,~,~]=scalar_newton_sec(U_o+1e-7*tang_init,D1,L,N,sec,tol,max_it,x,verbose,solution,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);






for j=1:length(parseq)
    
    if parseq(j)==1
        par=1;
mu1=[];
a1=[];
b1=[];
c1=[];
ds=1e-3;%initial step length
max_it=1000; % max number of Newton iterations
num=2000000; % max number of continuation steps 

break_sign=false;

for cont_steps=1:num
    cont_steps

    if cont_steps==num
        disp('number of contiuation steps too low, transition point not found')
        return 
    end

    
    if break_sign
        break
    end
    
    if ds<1e-10
        disp('step length is less than 1e-10')
        break
    end
    
    sec=U_n-U_o;
    sec=sec/norm(sec);
    U_g=U_n+(ds*sec);

    [U_sol,status,n_steps]=scalar_newton_sec(U_g,D1,L,N,sec,tol,max_it,x,verbose,true,    Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);

    switch status
    %status:
    %1: no_solution found or exceed maximum iteration
    %2: solution found within maximum iteration
    %3: initial guess is bad
        case 1 %solution not found
            ds=ds/4;
            disp(['steps=',num2str(cont_steps)])
            disp(['decreasing step size since no solution found: ds=' num2str(ds)])

        case 2 %solution founded!
             figure(2)
             
             
                 up=U_n(N+3);
                 um=U_n(N+4);
                 nu=U_n(N+6);

               WU=um*Chi_m+up*Chi_p+ Chi_p.*(U_n(N+1)*x+U_n(N+2)).*exp(nu*x)+ U_n(1:N);
             


              plot(x,WU,'b');
              title(['Solution Shape with a=' num2str(U_n(N+1))  '   mu1='  num2str(U_n(N+7)) ])
%              drawnow


            
            U_o=U_n;
            U_n=U_sol;
            mu=U_sol(N+7);
            a=U_sol(N+1);
            b=U_sol(N+2);
            c=U_sol(N+5);
 
            %fetch the critical mu
            if ~isempty(mu1) && (a)*(a1(end))<0
                break_sign=true;
               
            end
            
            c1=[c1;c];
            mu1=[mu1;mu];
            a1=[a1;a];
            b1=[b1;b];
            if n_steps>7%bad guess
                ds=ds/2;
                disp(['steps=',num2str(cont_steps)])
                disp(['decreasing step size since for doing multiple iteration: ds=' num2str(ds)])
            end 
            
            if n_steps<5 && ds<1e-2 
            ds=ds*1.5;
            disp(['steps=',num2str(cont_steps)])
            disp(['increasing step size for good initial guess: ds=' num2str(ds)])
            end

            
        otherwise %case 3
            ds=ds/2;
            disp(['steps=',num2str(cont_steps)])
            disp(['decreasing step size since bad initial guess: ds=' num2str(ds)])           
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% cont transition

    end
    
    if parseq(j)==2

U_good=U_n;


par=2;
num=1;
direction= 1;


a2=[];
b2=[];
mu2=[];
c2=[];
ds=1e-1;%initial step length
max_it=100;
break_sign=false;




for cont_steps=1:num
    
    if break_sign
        break
    end
    
    if ds<1e-10
        disp('step length is less than 1e-10')
        break
    end
    
    sec=U_n-U_o;
    sec=sec/norm(sec);
    U_g=U_n+(ds*sec);

    [U_sol,status,n_steps]=scalar_newton_sec(U_g,D1,L,N,sec,tol,max_it,x,verbose,true,    Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);

    switch status
    %status:
    %1: no_solution found or exceed maximum iteration
    %2: solution found within maximum iteration
    %3: initial guess is bad
        case 1 %solution not found
            ds=ds/4;
            disp(['steps=',num2str(cont_steps)])
            disp(['decreasing step size since no solution found: ds=' num2str(ds)])

        case 2 %solution founded!
             figure(2)
             
                 up=U_n(N+3);
                 um=U_n(N+4);
                 nu=U_n(N+6);

               WU=um*Chi_m+up*Chi_p+ Chi_p.*(U_n(N+1)*x+U_n(N+2)).*exp(nu*x)+ U_n(1:N);

             
             plot(x,WU,'b');
             title(['Solution Shape with a=' num2str(U_n(N+1))  '   mu='  num2str(U_n(N+7)) ])
             drawnow


            
            U_o=U_n;
            U_n=U_sol;
            mu=U_sol(N+7);
            a=U_sol(N+1);
            b=U_sol(N+2);
            c=U_sol(N+5);
                        

            

           
            
            mu2=[mu2;mu];
            a2=[a2;a];
            c2=[c2;c];
            b2=[b2;b];
            if n_steps>7%bad guess
                ds=ds/2;
                disp(['steps=',num2str(cont_steps)])
                disp(['decreasing step size since for doing multiple iteration: ds=' num2str(ds)])
            end 
            
            if n_steps<5 && ds<1e-2 
            ds=ds*1.5;
            disp(['steps=',num2str(cont_steps)])
            disp(['increasing step size for good initial guess: ds=' num2str(ds)])
            end

            
        otherwise %case 3
            ds=ds/2;
            disp(['steps=',num2str(cont_steps)])
            disp(['decreasing step size since bad initial guess: ds=' num2str(ds)])           
    end
end

    end
    

    
    if parseq(j)==4

            par=4;
mu4=[];
a4=[];
c4=[];
b4=[];
ds=1e-3;%initial step length
max_it=100;
num=200;
break_sign=false;

for cont_steps=1:num
    cont_steps
    
    if break_sign
        break
    end
    
    if ds<1e-10
        disp('step length is less than 1e-10')
        break
    end
    
    sec=U_n-U_o;
    sec=sec/norm(sec);
    U_g=U_n+(ds*sec);

    [U_sol,status,n_steps]=scalar_newton_sec(U_g,D1,L,N,sec,tol,max_it,x,verbose,true,    Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);

    switch status
    %status:
    %1: no_solution found or exceed maximum iteration
    %2: solution found within maximum iteration
    %3: initial guess is bad
        case 1 %solution not found
            ds=ds/4;
            disp(['steps=',num2str(cont_steps)])
            disp(['decreasing step size since no solution found: ds=' num2str(ds)])

        case 2 %solution founded!
             figure(2)
             
             
                 up=U_n(N+3);
                 um=U_n(N+4);
                 nu=U_n(N+6);

               WU=um*Chi_m+up*Chi_p+ Chi_p.*(U_n(N+1)*x+U_n(N+2)).*exp(nu*x)+ U_n(1:N);
             


             plot(x,WU,'b');
             title(['Solution Shape with c=' num2str(U_n(N+5))  '   mu='  num2str(U_n(N+7)) ])
             drawnow


            
            U_o=U_n;
            U_n=U_sol;
            mu=U_sol(N+7);
            a=U_sol(N+1);
            c=U_sol(N+5);
            b=U_sol(N+2);
 
           
            c4=[c4;c];
            mu4=[mu4;mu];
            a4=[a4;a];
            b4=[b4;b];
            if n_steps>7%bad guess
                ds=ds/2;
                disp(['steps=',num2str(cont_steps)])
                disp(['decreasing step size since for doing multiple iteration: ds=' num2str(ds)])
            end 
            
            if n_steps<5 && ds<1e-2 
            ds=ds*1.5;
            disp(['steps=',num2str(cont_steps)])
            disp(['increasing step size for good initial guess: ds=' num2str(ds)])
            end

            
        otherwise %case 3
            ds=ds/2;
            disp(['steps=',num2str(cont_steps)])
            disp(['decreasing step size since bad initial guess: ds=' num2str(ds)])           
    end
end


    end
end
        
        



