function [p e] = lm_opt(p0,fnc,grad,xdata,ydata)
    N = length(xdata);
    Np = length(p0);
    max_iters = 1e7;
    tol = 1e-4;
    lambda = 0.1;
    updateJ=1;
    J = zeros(N,Np);
    p = p0';
    
    for it = 1:max_iters
        if updateJ == 1
            %Jacobian
            for i=1:N
                J(i,:)=grad(p,xdata(i));
            end
            d = ydata - fnc(p,xdata);
            %Hessian approx
            H = transpose(J)*J;
            if it == 1
                e = dot(d,d);
            end
        end
        
        %Add dampening factor
        H = H+(lambda*eye(Np,Np));
        
        %Adjust arg vector
        dp = -inv(H)*(J'*d(:));
        p_tmp = p+dp;
        
        %Find error
        d_tmp = ydata - fnc(p_tmp,xdata);
        e_tmp = dot(d_tmp,d_tmp);
        %e-e_tmp
        
        %If error decreases, update parameters with tmp set and decrease
        %dampening factor. Otherwise try with bigger dampening factor
        if e_tmp < e
            lambda = lambda*0.1;
            p = p_tmp;
            e = e_tmp;
            %disp(e);
            updateJ=1;
        else
            updateJ = 0;
            lamda = lambda*10;
        end
        if (e < tol | norm(dp) < tol)
            disp('Exit from Convergence')
            break;
        end
    end
end