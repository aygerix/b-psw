import re, sys
src = open(sys.argv[1]).read()
m = re.search(r"ring\s+R\s*=\s*0\s*,\s*\(([^)]*)\)\s*,\s*dp\s*;", src)
vars_ = [v.strip() for v in m.group(1).split(",")]
body = re.search(r"ideal\s+I\s*=(.*?);", src, re.S).group(1)
gens = [g.strip() for g in body.strip().split(",\n")]
assert all(g for g in gens), "empty generator"
out = []
out.append("// k4q.sing translated to Calyx (Magma language) by sing2calyx.py. Run: calyx < k4q.m")
out.append("// Singular: ring R=0,(%s),dp; option(redSB); std(I); dim; size." % ",".join(vars_))
out.append('P<%s> := PolynomialRing(Rationals(), %d, "grevlex");' % (",".join(vars_), len(vars_)))
out.append("I := ideal<P |")
out.append(",\n".join(gens))
out.append(">;")
out.append("printf \"generators: %o, variables: %o\\n\", #Basis(I), Rank(P);")
out.append("t0 := Cputime(); r0 := Realtime();")
out.append("G := GroebnerBasis(I);")
out.append("printf \"GroebnerBasis: %o s cpu, %o s real\\n\", Cputime(t0), Realtime(r0);")
out.append("t1 := Cputime();")
out.append("d := Dimension(I);")
out.append("printf \"Dimension: %o s cpu\\n\", Cputime(t1);")
out.append('print "dim:"; print d; print "size:"; print #G;')
out.append('F := "%s";' % sys.argv[3])
out.append('PrintFile(F, "// reduced Groebner basis of k4q (grevlex, " cat Sprint(#G) cat " elements), Krull dimension " cat Sprint(d) : Overwrite := true);')
out.append("for g in G do PrintFile(F, g); end for;")
out.append("quit;")
open(sys.argv[2], "w").write("\n".join(out) + "\n")
print(len(vars_), "variables,", len(gens), "generators")
