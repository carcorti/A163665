\\\\ Independent PARI/GP verification of the consolidated A163665 terms.
\\\\ This script verifies all 28 positive hits, not the negative range.

hits = [ \
[1, 2], \
[2, 3], \
[8, 19], \
[12, 37], \
[128, 719], \
[240, 1511], \
[720, 5443], \
[6912, 69709], \
[32768, 386093], \
[142560, 1907819], \
[712800, 10777931], \
[1140480, 17819101], \
[1190400, 18653749], \
[3345408, 56125547], \
[3571200, 60163267], \
[5702400, 98911811], \
[14859936, 272887613], \
[29719872, 567611663], \
[50319360, 989060309], \
[118879488, 2444540149], \
[2147483648, 50685770167], \
[3889036800, 94206075193], \
[4389396480, 106882008719], \
[21946982400, 571264143487], \
[47416320000, 1272242879459], \
[92177326080, 2536961181761], \
[133145026560, 3715374607207], \
[331914240000, 9576704442589] \
];

ok = 1;

for (i = 1, #hits, \
    n = hits[i][1]; \
    k = hits[i][2]; \
    phi_sigma = eulerphi(sigma(n)); \
    prime_k = isprime(k); \
    row_ok = (phi_sigma == n && prime_k); \
    print(i, "\tn=", n, "\tk=", k, "\tphi_sigma_ok=", phi_sigma == n, "\tisprime_k=", prime_k, "\trow_ok=", row_ok); \
    if (!row_ok, ok = 0) \
);

print("PARI_GP_ALL_28_HIT_VERIFICATION=", if(ok, "OK", "FAIL"));

if (!ok, error("PARI/GP all-term hit verification failed"));

\\q
