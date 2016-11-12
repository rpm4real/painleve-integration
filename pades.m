function [a,b] = pades(z,u,up) %pade approximation centered at u(z) for step h
    %given initial conditions
    
    c(1) = up;
    
    % further coefficients found using mathematica, based on previous coefs
    c(2) = .5*(6*u^2+z);
    c(3) = (1/6)*(1 + 12*u*c(1));
    c(4) = (1/2)*(c(1)^2 + 2*u*c(2) );
    c(5) = (3/5)*(c(1)*c(2) + u*c(3) );
    c(6) = (1/5)*(c(2)^2 + 2*c(1)*c(3) + 2*u*c(4) );
    c(7) = (2/7)*(c(2)*c(3) + c(1)*c(4) + u*c(5) );
    c(8) = (3/28)*(c(3)^2 + 2*c(2)*c(4) + 2*c(1)*c(5) + 2*u*c(6) );
    c(9) = (1/6) *(c(3)*c(4) +c(2)* c(5)+c(1)* c(6)+u*c(7));
    c(10) = (1/15)* (c(4)^2+2 *c(3)* c(5)+2 *c(2) *c(6)+2* c(1)* c(7)+2 *u* c(8));
    c(11) = (6/55)* (c(4)* c(5) + c(3)* c(6) + c(2)* c(7) + c(1)* c(8) + u* c(9));
    c(12) = (1/22)* (2 *u *c(10)+c(5)^2+2 *c(4)*c(6)+2*c(3) *c(7)+2* c(2)*c(8)+2*c(1)* c(9));
    c(13) = (1/13)*(c(1) *c(10) + u* c(11)+ c(5) *c(6) + c(4)* c(7) + c(3)* c(8) + c(2)*c(9));
    c(14) = (3/91)* (2*c(1)*c(11)+2*u*c(12)+2*c(10)*c(2)+c(6)^2+2*c(5)*c(7)+2*c(4)*c(8)+2*c(3)*c(9));
    % we now have a formula for the 6-th order taylor expansion!
    
    nu = 7;
    
    %create the C matrix
    C = zeros(nu,nu);
    for i=0:nu-1
       for j=0:nu-1
           C(i+1,j+1) = c(nu+i-j);
       end
    end
    
    %create C vector
    Cvec = zeros(nu,1);
    for j = 1:nu
        Cvec(j) = c(nu+j);
    end
    Cvec = -Cvec;
    
    %Solve system for b coefficients
    b = C\Cvec;
    
    %second matrix of Cs
    C1 = diag(u*ones(nu,1),0);
    for d=1:nu-1    
       C1 = C1 + diag(c(d)*ones(nu-d,1),-d); 
    end
    
    %second C vector
    Cvec1 = zeros(nu,1);
    for i=1:nu
        Cvec1(i) = c(i);
    end
    
    % does matrix multiplication to find a coefficients
    a = C1*b + Cvec1;
    
    a = flipud(a);
    b = flipud(b);
    a = [u; a];
    b = [1;b];

end