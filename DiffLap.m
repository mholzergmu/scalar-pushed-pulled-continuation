function [D1,D2] = DiffLap(N,fourth_order,neumann,dx)

% Define Laplacian and differential operators
e = ones(N,1);


if fourth_order
    D2 = spdiags([(-1/12)*e (4/3)*e (-5/2)*e (4/3)*e (-1/12)*e], -2:2, N, N); %fourth order needs dx~.04, not smaller
    D1=  spdiags([(1/12)*e (-2/3)*e 0*e (2/3)*e (-1/12)*e], -2:2, N, N);
  

  if neumann
     


       cf=10/12;
       cc=-1/12;

    D2(1,1)=-15/12+cf*(-10/3);
    D2(1,2)=-4/12+cf*6;
    D2(1,3)=14/12-cf*2;
    D2(1,4)=-6/12+cf/3;
    D2(1,5)=1/12;
    
    D2(2,1)=D2(2,1)+cc*(-10/3);
    D2(2,2)=D2(2,2)+cc*6;
    D2(2,3)=D2(2,3)-cc*2;
    D2(2,4)=D2(2,4)+cc/3;
    
    D2(N,N)=-15/12+cf*(-10/3);
    D2(N,N-1)=-4/12+cf*6;
    D2(N,N-2)=14/12-cf*2;
    D2(N,N-3)=-6/12+cf/3;
    D2(N,N-4)=1/12;
    
    D2(N-1,N)=D2(N-1,N)+cc*(-10/3);
    D2(N-1,N-1)=D2(N-1,N-1)+cc*6;
    D2(N-1,N-2)=D2(N-1,N-2)-cc*2;
    D2(N-1,N-3)=D2(N-1,N-3)+cc/3;
    
    D1(1,2)=0;
    D1(1,1)=0;
    D1(1,3)=0;
    D1(N,N-3)=0;
    D1(N,N-2)=0;
    D1(N,N-1)=0;
    D1(N,N)=0;
    
    D1(2,1)=-3/12;
    D1(2,2)=-10/12;
    D1(2,3)=18/12;
    D1(2,4)=-6/12;
    D1(2,5)=1/12;
    D1(N-1,N)=-3/12;
    D1(N-1,N-1)=-10/12;
    D1(N-1,N-2)=18/12;
    D1(N-1,N-3)=-6/12;
    D1(N-1,N-4)=1/12;
   
  end
else 
    D2 = spdiags([e -2*e e], -1:1, N, N);
  D1=spdiags([-3/2*e 2*e -1/2*e ],0:2,N,N);
    if neumann
        D2(1,1)=-1;D2(N,N)=-1;
        D1(1,1)=-1/2; D1(N,N)=-1/2;
    end
end




D1 = D1/dx; D2 = D2 / dx^2; 
end

