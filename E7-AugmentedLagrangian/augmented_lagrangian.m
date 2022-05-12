function y = augmented_lagrangian(A,B,W,gamma,x)
%%AUGMENTED_LAGRANGIAN applies the augmented Lagrangian preconditioner
%%using direct solution of the inner blocks.

Agamma = A + gamma*B'*(W\B);
P = [Agamma B'; sparse(size(B,1),size(A,2)),-1/gamma*W];
y = P\x;

end