\\ quick local sanity run (tiny): Theorem U on each form
read("univ_lib.gp");
run("F",30,3,300,2000);
run("L",30,3,300,2000);
run("T",20,3,100,2000);
\\ classical controls
print("Chernick (6m+1)(12m+1)(18m+1), m=1: n=1729 sprp2=",sprp2(1729)," bpsw=",bpsw(1729));
print("composite F_p controls (Lemma 13a):"); forprime(q=3,200, my(F=fibonacci(q)); if(!ispseudoprime(F) && kronecker(5,F)==-1, print("  F_",q," slprp=",slprp(F)," sprp2=",sprp2(F)," D=",selD(F))));
print("composite Mersenne controls:"); forprime(q=3,80, my(M=2^q-1); if(!ispseudoprime(M), print1("  M_",q,":",sprp2(M)," ")));
print();
