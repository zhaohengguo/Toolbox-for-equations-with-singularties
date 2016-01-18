function [ Eigenvalue ] = cal_bessel( n )
% This function is used to calculate the eigenvalues and eigenfunctions
% of a Bessel-type equation. Equations of this kind always contain a
% singularity, thus using trigonometric function as the basis may not work.
% Instead, we use Bessel functions as the basis to expand this system and
% avoid these singularities.

% n: the power of the interaction potential x^n.
% M_size: the dimension of the basis used to expand the system.
% A: the size of the interval used to calculating the integrals. The value
% of a should be large enough that the value of every base function at a
% is small enough.
% lim: a given upper limit to keep A not too big. We need to bear
% in mind that the program may fail to work if the interval of integrals
% is too large..
% Eigenvalue_one_d : the unsorted eigenvalues of the equation.
% Eigenvalue: the eigenvalues after being sorted. Listed from the smallest
% to the largest.
% Eigenvector_in_order: the eigenvectors of the equation. The first one
% corresponds to the smallest eigenvalue, and the last one corresponds to
% the largest.

load roots_1.mat; % We need to know where the basis functions come to zero.

M_size = 100;

lim = 1000000;
A = lim^(1/(2*n-2));
% a=10; % We can directly set a=10 if n is not too large.

module = zeros(1, M_size);

for i = 1:M_size
    fun = @(x) x.*besselj(1,1./A*roots_1(1, i)*x).*besselj(1,1./A*roots_1(1, i)*x);
    module(1, i) = integral(fun, 0, A);  
end;


M = zeros(M_size, M_size);

for i = 1:M_size
    for j = 1:M_size
        fun = @(x) (x.*(1./A*roots_1(1,i))^2 + x.*((n-2)*x.^(n-2)+x.^(2*n-2))).*besselj(1,1./A*roots_1(1, i)*x).*besselj(1,1./A*roots_1(1, j)*x);
        M(i, j) = 1/sqrt(module(1,i)*module(1,j))*integral(fun, 0, A);
    end;
end;

[Eigenvector_one_d, Eigenvalue_one_d_m] = eig(M);

Eigenvalue_one_d = zeros([1, M_size]);

for i =1:M_size
    Eigenvalue_one_d(i) = Eigenvalue_one_d_m(i,i);
end;

Eigenvalue = sort(Eigenvalue_one_d);

Eigenvector_in_order = zeros(M_size);

for i=1:M_size
    [A, B] = find(Eigenvalue_one_d_m==Eigenvalue(i));
    Eigenvector_in_order(1:M_size, i) = Eigenvector_one_d(1:M_size, A);
end;


end

