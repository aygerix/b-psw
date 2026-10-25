\\ Theorem U at scale + positive controls (VM). Output: runU.out
read("univ_lib.gp");
default(parisizemax, 1000000000);
run("F",80,3,20000,20000);
run("F",60,4,3000,20000);
run("F",40,5,1000,20000);
run("F",30,6,300,20000);
run("L",80,3,20000,20000);
run("L",40,5,1000,20000);
run("T",60,3,3000,20000);
run("T",30,5,300,20000);
