import re, sys
s=open(sys.argv[1]).read(); ch=sys.argv[3] if len(sys.argv)>3 else '0'
vs=re.search(r"ring R=0,\((.*?)\),dp",s).group(1)
body=re.search(r"ideal I=\n(.*?);\n",s,re.S).group(1)
gens=[g.strip().replace('\n','') for g in body.split(',\n')]
open(sys.argv[2],'w').write(vs+"\n"+ch+"\n"+",\n".join(gens)+"\n")
print(len(vs.split(',')),"vars",len(gens),"polys")
