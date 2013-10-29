% Copyright, 2011, Kenny Erleben, DIKU.
clear all;
close all;
clc;

lambda1 = 1.0;
max_iter = 50;
tol_rel  = 0.0;
tol_abs  = 0.0;
num_variables = 8;


A = [    4.0006         1        -2         1       1.5       1.5       1.5       1.5
           1    4.0006         1        -2       1.5       1.5       1.5       1.5
          -2         1    4.0006         1      -1.5      -1.5      -1.5      -1.5
           1        -2         1    4.0006      -1.5      -1.5      -1.5      -1.5
         1.5       1.5      -1.5      -1.5    4.0006         1         1         4
         1.5       1.5      -1.5      -1.5         1    4.0006         4         1
         1.5       1.5      -1.5      -1.5         1         4    4.0006         1
         1.5       1.5      -1.5      -1.5         4         1         1    4.0006 ];
         
b1  =  [  -1
    -1
    -1
    -1
		0
		0
		0
		0 ];

b  =  [  -1
    -1
    -1
    -1
		-1
		-1
		-1
		-1 ];

x	 = zeros(num_variables,1);


x0     = zeros(num_variables,1);
    
%[z1 e1 i1 f1 conv1 m1] = psor(A, b, x0, lambda1, max_iter, tol_rel, tol_abs, true);
display('lemke no tangential force (friction should be zero)')
[z1 err] = lemke(A,b1,x0);
display(z1);

display('lemke with tangential force (non-zero friction)')
[z1 err] = lemke(A,b,x0);
display(z1);

%
lo = [ 0.
			0.
			0.
			0.
			-0.1
			-0.1
			-0.1
			-0.1 ];

hi = [ 10000.
			10000.
			10000.
			10000.
			0.1
			0.1
			0.1
			0.1 ];


