import mpmath as mp, numpy as np, json, pathlib
mp.mp.dps=55
z2,z3,z4=[mp.zeta(s) for s in (2,3,4)]
C=1/(z3-1);Z=(z3-1)/(z3-z4);c=4*Z
fpp=C/Z**2*(z2-3*z3+2*z4);D=1+Z**3*fpp/2
rows=[]
def I(j):
 if j==0:return mp.mpf(1)
 return Z*C*(-z3/(j*(j+1))+z2/j**2-(mp.digamma(j+1)+mp.euler)/j**3)
for j in [1,3,10,100,1000,10000]:
 row={'j':j,'scaled_integral_ratio':float(j*j*I(j)/(2*D))}
 if j<=10:
  quad=mp.quad(lambda r:r**j*Z*C*(z3-mp.polylog(3,r)/r) if r else mp.mpf(0),[0,mp.mpf('.8'),1])
  row['quadrature_abs_error']=float(abs(quad-I(j)))
 rows.append(row)
N=100000
# d recurrence gives exact central binomial normalization without large integers.
d=np.empty(N+202);d[0]=1.
for k in range(1,len(d)):d[k]=d[k-1]*(2*k-1)/(2*k)
h=2*np.arange(len(d))*d
ks=np.arange(2,N+2,dtype=float)
positive=float(C)/(2*ks**4*d[2:N+2])
# Harmonic/transition-row check: negative sum is finite at any fixed perimeter.
negative={j:float(2/c*d[j-1]*float(I(j-1))) for j in range(1,102)}
checks=[]
for p in [1,2,10,100]:
 pp=np.arange(1,N+1)
 plus=np.sum(positive*h[p+pp]/h[p])
 minus=sum(negative[j]*h[p-j]/h[p] for j in range(1,p))
 checks.append({'p':p,'row_sum_truncated':float(plus+minus),'tail_cutoff':N})
res={'C':float(C),'Z':float(Z),'c':float(c),'fpp':float(fpp),'D':float(D),'critical_residual':float(Z*C*(z3-z4)-1),'D_identity_residual':float(D-Z*C*(z2-z3)/2),'integral_checks':rows,'transition_rows':checks,'max_memory_MB_estimate':15}
pathlib.Path('experiment_results.json').write_text(json.dumps(res,indent=2),encoding='utf8')
print(json.dumps(res,indent=2))
