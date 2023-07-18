function F= scalar_sec( U,U_init,sec,D1,L,N,x,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);

w=U(1:N);
a=U(N+1); %a
b=U(N+2); %b
up=U(N+3);
um=U(N+4);
c=U(N+5);
nu=U(N+6);
mu=U(N+7);



    u0=um*Chi_m+up*Chi_p;

    %refined ansatz and derivatives
    
      
    
       
       
    re=(a*x+b).*exp(nu*x); % far-field term, a solution to the linear equation near the unstable state
    
    
    Wu=u0 + re.*Chi_p + w;  % the u component
    
    dxWu=D1*Wu;  % the derivative of the u component
    
    
    % 
    
    Fu=L*w+c*D1*w+commLChip*re+c*(D1_Chi_p).*re+f(Wu,dxWu,mu)+...
        +um*(L_Chi_m+c*D1_Chi_m)+up*(L_Chi_p+c*D1_Chi_p)-duf(0,0,mu)*Chi_p.*re;
    
    
      
    F1=sum(exp(-x.^2).*(Wu-(up+um)/2)); % phase condition 
    if par~=2 
    F2=sec'*(U-U_init); % secant continuation condition 
    elseif par==2 % a=0 located pushed to pulled transition 
        F2=a;
    end
    F3= w(end)+w(end-1);  % transverse boundary condition for core function 
    F4=f(up,0,mu); % steady state condition ahead of front
    F5=f(um,0,mu); % steady state condition behind front 
    F6=d(nu,c,mu); % dispersion relation near unstable state
    if par~=4
    F7=dn(nu,c,mu);  % unless pushed front continuation is taking place, impose double root condition 
    elseif par==4
        F7=a;
    end
    
    
    
    
 
    
    F=[Fu;F1;F2;F3;F4;F5;F6;F7];
 
    
end
