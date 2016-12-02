
%By using this function, you can solve larger systems more efficiently as there is no need to store the entire matrix A:

n2 = 21;
b2 = applyMoler(ones(n2,1));
tol = 1e-6;
maxit = 15;
M2 = spdiags((1:n2)',0,n2,n2);
[x2,flag2,rr2,iter2,rv2] = pcg(@applyMoler,b2,tol,maxit,M2);
% Using pcg with a Preconditioner
% This example demonstrates how to use a preconditioner matrix with pcg.
% 
% Create an input matrix and try to solve the system with pcg.

A = delsq(numgrid('S',100));
b = ones(size(A,1),1);
xi=zeros(size(A,1),1);
[x0,fl0,rr0,it0,rv0] = pcg(A,b,1e-8,100,[],[],xi);
% fl0 is 1 because pcg does not converge to the requested tolerance of 1e-8 within the requested maximum 100 iterations. A preconditioner can make the system converge more quickly.
% 
% Use ichol with only one input argument to construct an incomplete Cholesky factorization with zero fill.

L = ichol(A);
[x1,fl1,rr1,it1,rv1] = pcg(A,b,1e-8,100,L,L',xi);
% fl1 is 0 because pcg drives the relative residual to 9.8e-09 (the value of rr1) which is less than the requested tolerance of 1e-8 at the seventy-seventh iteration (the value of it1) when preconditioned by the zero-fill incomplete Cholesky factorization. rv1(1) = norm(b) and rv1(78) = norm(b-A*x1).
% 
% The previous matrix represents the discretization of the Laplacian on a 100x100 grid with Dirichlet boundary conditions. This means that a modified incomplete Cholesky preconditioner might perform even better.
% 
% Use the michol option to create a modified incomplete Cholesky preconditioner.

L = ichol(A,struct('michol','on'));
[x2,fl2,rr2,it2,rv2] = pcg(A,b,1e-8,100,L,L',xi);
% In this case you attain convergence in only forty-seven iterations.
% 
% You can see how the preconditioners affect the rate of convergence of pcg by plotting each of the residual histories starting from the initial estimate (iterate number 0).
[x11,h,hc01,hc02]=ICCG_0(A,b,xi,100,1e-8,L,0);

figure(500);
h1=semilogy(0:it0,rv0/norm(b),'b.');
hold on;
h2=semilogy(0:it1,rv1/norm(b),'r.');
h3=semilogy(0:it2,rv2/norm(b),'k.');
hold on
legend([h1,h2,h3,h],'No Preconditioner','IC(0)','MIC(0)','ICCG_{own}');
xlabel('iteration number');
ylabel('relative residual');
hold off;
