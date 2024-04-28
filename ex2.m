close all
clear all
clc

syms x y theta alfa beta gamma r l real;
q=[x y theta alfa beta gamma]';

A=[(sqrt(3)/2)*cos(theta)-0.5*sin(theta)  sin(theta)   (-sqrt(3)/2)*cos(theta)-0.5*sin(theta)
   0.5*cos(theta)+(sqrt(3)/2)*sin(theta)  -cos(theta)   0.5*cos(theta)-(sqrt(3)/2)*sin(theta)
                    l                        l                            l 
                    r                        0                            0
                    0                        r                            0
                    0                        0                            r                  ];

G=null(A');

g1=G(:,1);
g2=G(:,2);
g3=G(:,3);
G=simplify(G);

%First Lie Brackets
g4=jacobian(g2,q)*g1-jacobian(g1,q)*g2;     %[g1,g2]
g5=jacobian(g3,q)*g2-jacobian(g2,q)*g3;     %[g2,g3]
g6=jacobian(g3,q)*g1-jacobian(g1,q)*g3;     %[g1,g3]

disp('Accessibility rank condition:')
F=[G g4]
disp('Rank of F is:')
disp(rank(F))      %rank is increased to 4

fprintf('\n');
disp('Accessibility rank condition:')
F=[G g4 g5]
disp('Rank of F is:')
disp(rank(F))      %rank is increased to 5

fprintf('\n');
disp('Accessibility rank condition:')
F=[G g4 g5 g6]
disp('Rank of F is:')
disp(rank(F))     %rank is still 5 -> g6 is linearly dependent

%It's possible to explore all the combination of independent vectors of F,
%in order to find the last independent vector. The possible combinations
%are:
g7=jacobian(g4,q)*g1-jacobian(g1,q)*g4;     %[g1,g4]
g8=jacobian(g5,q)*g1-jacobian(g1,q)*g5;     %[g1,g5]
g9=jacobian(g4,q)*g2-jacobian(g2,q)*g4;     %[g2,g4]
g10=jacobian(g5,q)*g2-jacobian(g2,q)*g5;    %[g2,g5]
g11=jacobian(g4,q)*g3-jacobian(g3,q)*g4;    %[g3,g4]
g12=jacobian(g5,q)*g3-jacobian(g3,q)*g5;    %[g3,g5]
g13=jacobian(g5,q)*g4-jacobian(g4,q)*g5;    %[g4,g5]

fprintf('\n');
disp('Accessibility rank condition:')
F=[G g4 g5 g7 g8 g9 g10 g11 g12 g13];
disp(rank(F))  %rank is still 5 -> all the vectors are linearly depedent

fprintf('Maximum rank of F: %d\n', rank(F));      %the dimension of the accessibility distribution is 5