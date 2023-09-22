%% PHY 329 Homework 2
%% Ryan Schlimme (eid: rjs4499)
%% Problem 9.16

% Design function to solve pentadiagonal system without pivoting
% Test Case
A = [8 -2 -1 0 0;
    -2 9 -4 -1 0;
    -1 -3 7 -1 -2;
    0 -4 -2 12 -5;
    0 0 -7 -3 -15]
b = [5; 2; 1; 1; 5]


x = pentsolve(A,b)
A*x
%% Problem 10.10

A = [3 -2 1;
    2 6 -4;
    -1 -2 5];
b = [-10; 44; -26]

% The lu function automatically uses partial pivoting to avoid division by
% 0. Therefore, we can solve the system as follows.

[L,U] = lu(A);      % LU factorization with partial pivoting
d = L\b;            % intermediary solution
x = U\d             % actual solution

% Check
A*x
%% 
% Therefore, x1 = 1, x2 = 5, x3 = -3. We also checked that Ax = b and see our 
% result is correct.
%% Problem 10.11
%% Part c)

% Complete LU factorization without pivoting and compute the determinant
% using the result.

A = [8 2 1;
    3 7 2;
    2 3 9]

[L, U] = lu_nopivot(A)

L*U % Check factorization

detLU = det(U)*det(L)
det(A) % Check determinant
%% 
% I can confirm that the LU factorization without pivoting satisfied LU = A. 
% Additionally, by computing the product of the trace of the upper and lower diagonal 
% matrix, one finds that det(A) = 405.
%% Problem 10.12
% Given the following LU decomposition and vector:

L = [1 0 0;
    0.6667 1 0;
    -0.3333 -0.3636 1];

U = [3 -2 1;
    0 7.3333 -4.667;
    0 0 3.6364];

A = L*U;

b = [-10; 44; -26];
%% Part a)
% Compute the determinant

det(L)*det(U)
det(A) % check
%% 
% The determinant via LU decomposition is 80.0004.
%% Part b)
% Solve Ax = b

d = L\b;            % intermediary solution
x = U\d             % actual solution

A*x                 % check
%% 
% x = [0.999, 4.9997, -3.0004]
%% Problem 11.2
% Determine the matrix inverse for the following system

A = [-8 1 -2;
    2 -6 -1;
    -3 -1 7];

Ainv = inv(A)
%% 
% Our matrix inverse is given by Ainv.
%% Problem 11.6
% Determine the Frobenius (F), Column-Sum (1), and Row-Sum (infinity) norm of 
% A after scaling the matrix such that the maxiomum element in each row is equal 
% to 1.

A = [8 2 -10;
    -9 1 3;
    15 -1 6;]

A = [-0.8 -0.2 1;
    1 -1/9 -1/3;
    1 -1/15 6/15];

normF = norm(A, "fro")
norm1 = norm(A, 1)
normInf = norm(A, "inf")
%% 
% Frobenius Norm: 1.9920, Column-Sum Norm: 2.8000, Row-Sum Norm: 2
%% Problem 11.8
% Determine the spectral condition number (p = 2) for the system without normalization. 
% Compute condition number based on row-sum norm.

A = [1 4 9 16 25;
    4 9 16 25 36;
    9 16 25 36 49;
    16 25 36 49 64;
    25 36 48 64 81]

spec = cond(A, 2)
rowSum = cond(A, 1)
%% 
% The spectral condition number for the system is 8.037e16 and the row-sum condition 
% number is 1.723e18.
%% Problem 12.2
%% Part a)
% Use the Gauss-Seidel method to solve the system until the percent relative 
% error falls below 5%.

A = [0.8 -0.4 0;
    -0.4 0.8 -0.4;
    0 -0.4 0.8];

b = [41; 25; 105];

er = 0.05;

x = GaussSeidel(A, b, er)
%% 
% Using Gauss-Seidel with 5% error, we receive x = 167.8711, 239.1211, 250.8105.
%% Part b)
% Now use overrelaxation with lambda = 1.2

A = [0.8 -0.4 0;
    -0.4 0.8 -0.4;
    0 -0.4 0.8];

b = [41; 25; 105];

er = 0.05;

lambda = 1.2;

x = GaussSeidel(A, b, er, lambda)
%% 
% Using overrelaxation with lambda = 1.2, our results change slightly to 79.6875, 
% 150.9375, and 206.7188.
%% Problem 12.8
% Determine the solution to the nonlinear system using Newton-Raphson with initial 
% guesses x=y=1.2.

syms x;
syms y;

fval1 = -x^2+x+0.75-y;
fval2 = x^2-5*x*y-y;

fval3 = [fval1;
    fval2];

jac = [diff(fval1, x), diff(fval1, y);
    diff(fval2, x), diff(fval2, y)];

p1 = inv(jac)*fval3

syms tol;    %% Define the tolerance vector
syms err;    %% Define the error matix that checks how closely is the result of an iteration converging on the solution.
X02 = [1.2;  %% Define guess values for the variables
     1.2];
iter=10;       %% Define Number of Iterations. More the number of iterations, longer will the code take to converge on a solution.    
tol=[1e-6;
    1e-6];     %% Specify desired tolerances.
X2=X02;        
Xold2=X02;
for i=1:iter
    h12=subs(p1(1),{x,y},{X2(1),X2(2)}); %% Evaluation of the partial differentials with the guess values.
    h22=subs(p1(2),{x,y},{X2(1),X2(2)});
    u=[h12;
       h22];
    X2=X2-u;
    double(X2)
    err=abs(X2-Xold2);
    Xold2=X2;
    if (err(1)<tol(1) && err(2)<tol(2))
        break
    end
end

x = 1.3721
y = 0.2395

fval1 = -x^2+x+0.75-y
fval2 = x^2-5*x*y-y
%% 
% Using multivariate Newton-Raphson with 10 iterations maximum, we converge 
% to x = 1.3721 and y = 0.2395. If we check our answers, we find that fval1 and 
% fval2 are both approximately 0 which confirms our solution.
%% Problem 13.4
% After deriving the system of differential equations for a three mass, four 
% spring system as illustrated in the homework document, we can solve for the 
% eigenvalues and natural frequencies for the following values of spring constant 
% and mass. 

k1 = 15;
k4 = k1;
k2 = 35;
k3 = k2;

m = 1.5;

% Takes the form Ax = cx, x: eigenvector, c: eigenvalue
A = [-3.3333 2.3333 0;
    2.3333 -4.6667 2.3333;
    0 2.3333 -3.3333]

[u,v] = eig(A)

sqrt(7.3665)
sqrt(3.3333)
sqrt(0.6335)
%% Problem 13.5
% Since all values of mass are equal to 1, we can neglect the contribution of 
% m. Since w^2 is our eigenvalue, we remove that from the A matrix since the eig 
% function will subtract a constant (the eigenvalue) from the diagonal terms automatically.

k = 2;

A = [2*k -k -k;
    -k 2*k -k;
    -k -k 2*k];

[u,v] = eig(A)
sqrt(v(2,2))
%% 
% We find eigenvalues w^2 = 0, 6, 6. Therefore, our frequencies are 0 and 2.4495 
% rad/s with multiplicity 2.
%% Problem 13.6
% Since L = 1, our matrix simplifies.

C = 0.25

A = [1/C -1/C 0;
    -1/C 2/C -1/C;
    0 -1/C 2/C];

[u,v] = eig(A)
sqrt(v(1,1))
sqrt(v(2,2))
sqrt(v(3,3))
%% 
% By solving for the eigenvalues and taking the square root, we get resonant 
% frequencies of 0.8901, 2.4940, and 3.6039 rad/s. 
%% Define Our Functions

% LU factorization without pivoting (from stack overflow)
%%
function [L, U] = lu_nopivot(A)
n = size(A,1);
L = eye(n);
for k = 1:n
    L(k+1 : n, k) = A(k+1 : n,k) / A(k,k);
    for l = k+1:n
        A(l, :) = A(l, :) - L(l, k)*A(k,:);
    end
end
U = A;
end

% Pentadiagonal Solver
function x=pentsolve(A,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pentsolve.m
%
% Solve a pentadiagonal system Ax=b where A is a strongly nonsingular matrix
% 
% If A is not a pentadiagonal matrix, results will be wrong
%
% Reference: G. Engeln-Muellges, F. Uhlig, "Numerical Algorithms with C"
%               Chapter 4. Springer-Verlag Berlin (1996)
%
% Written by Greg von Winckel 3/15/04
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,N]=size(A);
% Check dimensions
if M~=N
    error('Matrix must be square');
    return; 
end
if length(b)~=M
    error('Matrix and vector must have the same number of rows');
    return;
end
x=zeros(N,1);
    
% Check for symmetry
if A==A'    % Symmetric Matrix Scheme
    
    % Extract bands
    d=diag(A);
    f=diag(A,1);
    e=diag(A,2);
        
    alpha=zeros(N,1);
    gamma=zeros(N-1,1);
    delta=zeros(N-2,1);
    c=zeros(N,1);
    z=zeros(N,1);
    
    % Factor A=LDL'
    alpha(1)=d(1);
    gamma(1)=f(1)/alpha(1);
    delta(1)=e(1)/alpha(1);
    
    alpha(2)=d(2)-f(1)*gamma(1);
    gamma(2)=(f(2)-e(1)*gamma(1))/alpha(2);
    delta(2)=e(2)/alpha(2);
    
    for k=3:N-2
        alpha(k)=d(k)-e(k-2)*delta(k-2)-alpha(k-1)*gamma(k-1)^2;
        gamma(k)=(f(k)-e(k-1)*gamma(k-1))/alpha(k);
        delta(k)=e(k)/alpha(k);
    end
    
    alpha(N-1)=d(N-1)-e(N-3)*delta(N-3)-alpha(N-2)*gamma(N-2)^2;
    gamma(N-1)=(f(N-1)-e(N-2)*gamma(N-2))/alpha(N-1);
    alpha(N)=d(N)-e(N-2)*delta(N-2)-alpha(N-1)*gamma(N-1)^2;
    
    % Update Lx=b, Dc=z
    
    z(1)=b(1);
    z(2)=b(2)-gamma(1)*z(1);
    
    for k=3:N
        z(k)=b(k)-gamma(k-1)*z(k-1)-delta(k-2)*z(k-2);
    end
    
    c=z./alpha;
        
    % Backsubstitution L'x=c
    x(N)=c(N);
    x(N-1)=c(N-1)-gamma(N-1)*x(N);
    
    for k=N-2:-1:1
        x(k)=c(k)-gamma(k)*x(k+1)-delta(k)*x(k+2);
    end
    
else        % Non-symmetric Matrix Scheme
    
    % Extract bands
    d=diag(A);
    e=diag(A,1);
    f=diag(A,2);
    h=[0;diag(A,-1)];
    g=[0;0;diag(A,-2)];
        
    alpha=zeros(N,1);
    gam=zeros(N-1,1);
    delta=zeros(N-2,1);
    bet=zeros(N,1);
    
    c=zeros(N,1);
    z=zeros(N,1);
           
    % Factor A=LR
    alpha(1)=d(1);
    gam(1)=e(1)/alpha(1);
    delta(1)=f(1)/alpha(1);
    bet(2)=h(2);
    alpha(2)=d(2)-bet(2)*gam(1);
    gam(2)=( e(2)-bet(2)*delta(1) )/alpha(2);
    delta(2)=f(2)/alpha(2);
    
    for k=3:N-2
        bet(k)=h(k)-g(k)*gam(k-2);
        alpha(k)=d(k)-g(k)*delta(k-2)-bet(k)*gam(k-1);
        gam(k)=( e(k)-bet(k)*delta(k-1) )/alpha(k);
        delta(k)=f(k)/alpha(k);
    end
    
    bet(N-1)=h(N-1)-g(N-1)*gam(N-3);
    alpha(N-1)=d(N-1)-g(N-1)*delta(N-3)-bet(N-1)*gam(N-2);
    gam(N-1)=( e(N-1)-bet(N-1)*delta(N-2) )/alpha(N-1);
    bet(N)=h(N)-g(N)*gam(N-2);
    alpha(N)=d(N)-g(N)*delta(N-2)-bet(N)*gam(N-1);
    % Update b=Lc
    c(1)=b(1)/alpha(1);
    c(2)=(b(2)-bet(2)*c(1))/alpha(2);
    
    for k=3:N
        c(k)=( b(k)-g(k)*c(k-2)-bet(k)*c(k-1) )/alpha(k);
    end
    
    
    % Back substitution Rx=c
    x(N)=c(N);
    x(N-1)=c(N-1)-gam(N-1)*x(N);
    for k=N-2:-1:1
        x(k)=c(k)-gam(k)*x(k+1)-delta(k)*x(k+2);    
    end
   
end
end

% Gauss-Seidel Method

function x = GaussSeidel(A, b, es, maxit, lambda)
if nargin<2, error('at least 2 input arguments required'), end
if nargin<4 || isempty(maxit), maxit=50;end
if nargin<3 || isempty(es), es = 0.00001; end
if nargin<5 || isempty(relax), lambda = 1; end
[m,n] = size(A);
if m~=n, error('Matrix A must be square');end
C = A;
for i=1:n
    C(i,i) = 0
    x(i) = 0;
end
x = x';
for i=1:n
    C(i, 1:n) = C(i, 1:n)/A(i,i);
end
for i=1:n
    d(i) = b(i)/A(i,i);
end
iter = 0;
while (1)
    xold = x;
    for i = 1:n
        x(i) = d(i) - C(i, :)*x;
        x(i) = lambda*x(i) + (1- lambda)*xold(i);
        if x(i) ~= 0
            ea(i) = abs((x(i) - xold(i))/x(i));
        end
    end
    iter = iter + 1;
    if max(ea)<=es || iter>= maxit, break, end
end
end


% Multivariate Newton-Raphson
function [x, f, ea, iter] = newtmult(func, x0, es, maxit, varargin)
% func: function vector, x0: initial guess, es: error, maxit:
% maximum iterations
if nargin<2, error('at least 2 input args requires'), end
if nargin<3 || isempty(es), es = 0.0001; end
if nargin<4 || isempty(maxit), maxit = 50; end
iter = 0;
x = x0;
while(1)
    [J,f] = func(x, varargin{:});
    dx = J\f;
    x = x-dx;
    iter = iter + 1;
    ea = 100*max(abs(dx./x));
    if iter>= maxit || ea<=es, break, end
end
end