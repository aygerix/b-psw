{c=0; bad=0; first=List();
for(Y=3,19999, a=Y^2-Y-1; b=Y^3-Y^2-2*Y+1; d=Y^3-3*Y-1;
  if(isprime(a)&&isprime(b)&&isprime(d), c++; n=a*b*d;
    if((n-1)%(a-1) || (n-1)%(b-1) || (n-1)%(d-1), bad++);
    if(#first<3, listput(first,[Y,a,b,d]))));
print("all-prime Y in [3,20000): ",c,"  Korselt failures: ",bad);
print(first);}
