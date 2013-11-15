% A0,b0 is the LCP with no tangential force, so friction should be zero
% A1,b1 is the LCP with some tangential force, so friction should be non-zero
% The idea is to introduce lower and upper bound for the friction
% and convert those bounds into the A and b, so that friction is clamped
% Appendix A1 of this paper has a derivation, see http://www.cs.duke.edu/~parr/nips10.pdf

clear all;
close all;
clc;

lambda1 = 1.0;
max_iter = 10;
tol_rel  = 0.0;
tol_abs  = 0.0;
num_variables = 8;

A0 =   [    4.0006         1        -2         1       1.5       1.5       1.5       1.5
           1    4.0006         1        -2       1.5       1.5       1.5       1.5
          -2         1    4.0006         1      -1.5      -1.5      -1.5      -1.5
           1        -2         1    4.0006      -1.5      -1.5      -1.5      -1.5
         1.5       1.5      -1.5      -1.5    4.0006         1         1         4
         1.5       1.5      -1.5      -1.5         1    4.0006         4         1
         1.5       1.5      -1.5      -1.5         1         4    4.0006         1
         1.5       1.5      -1.5      -1.5         4         1         1    4.0006 ];
     
b0  =  [  -1
    -1
    -1
    -1
		0
		0
		0
		0 ];

A1 =     [    4.0006         1        -2         1       1.5       1.5       1.5       1.5
           1    4.0006         1        -2      -1.5      -1.5      -1.5      -1.5
          -2         1    4.0006         1      -1.5      -1.5      -1.5      -1.5
           1        -2         1    4.0006       1.5       1.5       1.5       1.5
         1.5      -1.5      -1.5       1.5    4.0006         4         1         1
         1.5      -1.5      -1.5       1.5         4    4.0006         1         1
         1.5      -1.5      -1.5       1.5         1         1    4.0006         4
         1.5      -1.5      -1.5       1.5         1         1         4    4.0006 ];
     
b1  =  [  -1
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
%display('lemke no tangential force (friction should be zero)')
%[z1 err] = lemke(A0,b0,x0);
%display(z1);

%display('lemke with tangential force (non-zero friction)')
%[z1 err] = lemke(A1,b1,x0);
%display(z1);

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

%Method 0, not working yet



%--- Definition of bounded LCP
%
% y = A x + b
%
% x = l     => y >= 0
% x = h     => y <= 0
% l < x < h => y  = 0
%
% Split y into positive and negative components y = y^+ - y^-
%
% Then we rewrite as follows
%
%    y^+ > 0 => x-l = 0
%    y^+ = 0 => x-l > 0
%    y^- > 0 => h-x = 0
%    y^- = 0 => h-x > 0
%
% From y = A x + b and B*A = I we have
%
%    y - b = A x
%        x = -B b + B y
%        x = -B b + B y^+ - B y^-
%
% Substitution yields
%
% y^+ > 0 =>   B y^+  - B y^- + (-B*b-l) = 0
% y^+ = 0 =>   B y^+  - B y^- + (-B*b-l) > 0
% y^- > 0 => - B y^+  + B y^- + (h+B*b) = 0
% y^- = 0 => - B y^+  + B y^- + (h+B*b) > 0
%
% Introducing matrix notation we may write
%
% |x^+ | =   |  B  -B|  | y^+ |   | (-B*b - l) |
% |x^- | =   | -B   B|  | y^- | + | (h + B*b) |
%    w   =       M         z    +       q
%
% Then we simply have  0<= w  compl. z >= 0
%

B = pinv(A1);
M = [B  -B;
    -B   B];
q = [ (-B*b1 - lo)'   (hi + B*b1)' ]';
z0 = zeros(16,1);
%[z1 err] = lemke(M, q, z0);
%[z1 e1 i1 f1 conv1 m1] = psor(M, q, z0, lambda1, max_iter, tol_rel, tol_abs, true);
[z1 e1 i1 f1 conv1 m1] = fischer_newton(M, q, z0, max_iter, tol_rel, tol_abs, 'random', true);

w1       = M * z1 + q;
display('Method 0: lemke with tangential force, friction clamped to lo=-0.1 hi=0.1]')
display(z1);
display(w1);
display('Complementarity test of LCP')
display(w1'*z1);

%--- Converting back from LCP solution space to BLCP solution space
y_plus  = z1(1:8);
y_minus = z1(9:16);
y1      = y_plus - y_minus;
x1      = B*(y1-b1);

display('Solution of BLCP')
display(x1);
display(y1);

display('Complementarity test of BLCP')
H = min(hi-x1,min(x1-lo,-y1));
display(H'*H);

%
%Method 1, convert the problem from a LCP to a QP and added bound constraints 
%on the variables, which yields the solution. 
%Effectively, if you wish to produce the "more physical" result, 
%then you would only constrain the normal directions of A:
display ('quadprog')
quadprog(A1,b1*0.5,-A1(1:4,:), b1(1:4), [], [], lo, hi)



