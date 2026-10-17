\\ symmetric unit pairs: Res(Phi_a+1,Phi_b)=+-1 AND Res(Phi_b+1,Phi_a)=+-1
N=160; P=vector(N,d,polcyclo(d)); U=matrix(N,N);
for(a=1,N, A=P[a]+1; for(b=1,N, if(b!=a, U[a,b]=(abs(polresultant(A,P[b]))==1))));
s=0; for(a=1,N, for(b=a+1,N, if(U[a,b]&&U[b,a], s++; print("sym pair ",a," ",b))));
print("symmetric pairs: ",s);
\\ Lucas analogue: Psi_d = min poly of 2cos(2pi/d); need Res(Psi_a-1,Psi_b)=+-1 both ways
Ps=vector(N,d, if(d<=2, x-(d==1)*2+(d==2)*2, my(f=polresultant(polcyclo(d,y), y^2-x*y+1, y)); f=factor(f)[1,1]; f/pollead(f)));
s=0; for(a=1,N, for(b=a+1,N, if(abs(polresultant(Ps[a]-1,Ps[b]))==1 && abs(polresultant(Ps[b]-1,Ps[a]))==1, s++; print("Lucas sym pair ",a," ",b))));
print("Lucas symmetric pairs: ",s);
