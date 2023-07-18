function J= scalar_dF_sec(U,U_init,sec,D1,L,N,x,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);


u=U(1:N);
a=U(N+1); 
b=U(N+2); 
up=U(N+3);
um=U(N+4);
c=U(N+5);
nu=U(N+6);
mu=U(N+7);



    u0=um*Chi_m+up*Chi_p;

    %refined ansatz and derivatives
    
    
    re=(a*x+b).*exp(nu*x);
    
    
    Wu=u0 + re.*Chi_p + u;  
    dxWu=D1*Wu;

    


    eps=1e-4;
    


 
    
    
    Juu=L+c*D1+spdiags(duf(Wu,dxWu,mu),0,N,N)+spdiags(duxf(Wu,dxWu,mu),0,N,N)*D1;    
    

    Ja=Juu;  % NxN
    Ja=[Ja;exp(-x.^2)']; % F1
    if par~=2
    Ja=[Ja;sec(1:N)']; % F2
    elseif par==2
     Ja=[Ja;sparse(1,N)]; % F2
    end
    
    Ja(N+3,N)=1; % F3
    Ja(N+3,N-1)=1; %F3
    Ja=[Ja;sparse(4,N)]; % F4-F7

     
    

    
    


    F=@(U) scalar_sec(U,U_init,sec,D1,L,N,x,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);


    

    J1=(F(U+eps*[zeros(N,1);1;zeros(6,1)])-F(U))/eps;
    J2=(F(U+eps*[zeros(N+1,1);1;zeros(5,1)])-F(U))/eps;
    J3=(F(U+eps*[zeros(N+2,1);1;zeros(4,1)])-F(U))/eps;
    J4=(F(U+eps*[zeros(N+3,1);1;zeros(3,1)])-F(U))/eps;
    J5=(F(U+eps*[zeros(N+4,1);1;zeros(2,1)])-F(U))/eps;
    J6=(F(U+eps*[zeros(N+5,1);1;zeros(1,1)])-F(U))/eps;
    J7=(F(U+eps*[zeros(N+6,1);1])-F(U))/eps;  
   
    J=[Ja,J1,J2,J3,J4,J5,J6,J7];
    
    
end
