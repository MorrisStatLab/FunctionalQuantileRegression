function A_inv=invChol(A)

R=chol(A);

R_inv=inv(R);

A_inv=R_inv*R_inv';