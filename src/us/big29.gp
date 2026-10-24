mm=40666413944507969; OO=oddp(mm-1); print("m=",mm," factor m=",factor(mm)," O=",factor(OO));
t0=getabstime(); r=testm(29, mm, 2, 1); print("testm: ", r, "  time ", (getabstime()-t0)\1000, " s");
