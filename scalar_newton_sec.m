%status:
%1: no_solution found or exceed maximum iteration
%2: solution found within maximum iteration
%3: res_too_large: initial guess is bad
function [root,status,steps] = ...       
    scalar_newton_sec(U0,D1,L,N,sec,tol,max_it,x,verbose,solution,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn)
    

    steps = 1;
    root=U0;
    U_init=U0;      
    res=scalar_sec(U0,U_init,sec,D1,L,N,x,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);    

    if solution
        if norm(res)>1
            norm(res)
            status=3;
            return
        end
    end
%      
    while (norm(res,'inf') > tol) && (steps < max_it)
        J   = scalar_dF_sec(U0,U_init,sec,D1,L,N,x,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);
        steps = steps + 1;
        dU = -J\res;

        U0 = U0 + dU;
        res=scalar_sec(U0,U_init,sec,D1,L,N,x,Chi_p,D1_Chi_p,L_Chi_p,commLChip,Chi_m,D1_Chi_m,L_Chi_m,par,f,duf,duxf,d,dn);
        norm(res,'inf')
        if verbose
            residue_after_successful_newton=norm(res,'inf')
        end
    end
    
    if  norm(res,'inf') > tol
            status=1;
    else
        root = U0;
        status=2;
    end
    
    return
end
%  
