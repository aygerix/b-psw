\\ L_(2^k) (k >= 5, all composite) is AUTOMATICALLY a strong Lucas probable prime with Selfridge D = 5:
\\ L_(2^k) == 2 mod 5, every prime factor has z = 2^(k+1), and 2^(k+1) | L_(2^k) + 1 (induction:
\\ L_(2^(k+1)) + 1 = (L_(2^k) - 1)(L_(2^k) + 1)). So BPSW on L_(2^k) reduces to the strong base-2 test.
\\ This is the Lucas mirror of "every composite F_m is spsp(2)". Lottery: test k = 5..20 directly.
default(parisizemax, 4*10^9);
read("../fam/bpsw.gp");
{for(k=5,20, my(t=getabstime(), n=fibonacci(2^k-1)+fibonacci(2^k+1), D=selD(n), sl=slprp(n), sf=sprp2(n));
   printf("k=%d digits=%d D=%d strongLucas=%d spsp2=%d time %ds\n", k, #digits(n), D, sl, sf, (getabstime()-t)\1000));}
