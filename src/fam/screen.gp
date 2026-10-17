\\ Family F: n = Psi10(Y)*Psi14(Y)*Psi18(Y), with p_j - 1 | n - 1 identically in Y.
\\ Screen Y = L_M: Jacobi(5/n), weak Fermat base 2, weak Lucas U_{n+1}=0 (P=1,Q=-1).
P10(Y)=Y^2-Y-1;
P14(Y)=Y^3-Y^2-2*Y+1;
P18(Y)=Y^3-3*Y-1;
\\ sanity: identity check
{my(n=P10(y)*P14(y)*P18(y)); print("identity remainders: ", [(n-1)%(P10(y)-1), (n-1)%(P14(y)-1), (n-1)%(P18(y)-1)]);}
lucU(k,n)={my(M=Mod([1,1;1,0],n)^k); lift(M[1,2])};
{for(M=1,600, my(Y=fibonacci(M-1)+fibonacci(M+1), n=P10(Y)*P14(Y)*P18(Y), j=kronecker(5,n), f=(Mod(2,n)^(n-1)==1), l=0);
  if(j==-1, l=(lucU(n+1,n)==0));
  if(f||l, print("M=",M," jac=",j," fermat=",f," lucas=",l)));}
