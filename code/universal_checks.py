import numpy as np, json, pathlib, math
N=1050000
x=np.arange(1,N+1,dtype=float)
d=np.empty(N+1);d[0]=1;d[1:]=np.cumprod((2*x-1)/(2*x))
h=2*np.arange(N+1)*d
ms=np.array([1,4,64,4096,1000000],dtype=int)
models=[('atom_'+str(m),np.array([m]),np.array([1.])) for m in ms]
models += [('lacunary_equal',ms,np.ones(5)/5),('lacunary_skew',ms,np.array([.01,.01,.01,.01,.96]))]
res=[]
for name,m,weights in models:
 pos=weights/h[m+1];nu0=1-np.sum(weights)
 w=m*(2*m+1)*d[m]*pos
 for p in [1,2,4,8,16,64,256,1024,4096,16384]:
  k=np.arange(1,p,dtype=int)
  T=lambda a:np.sum(w/(m+a))
  rate=T(p)/math.sqrt(p)
  neg=d[k]/(2*k-1)*np.sum(w[None,:]/(m[None,:]+k[:,None]),axis=1) if p>1 else np.array([])
  probsneg=neg*h[p-k]/h[p]
  probspos=pos*h[p+m]/h[p]
  ys=np.concatenate((p-k,p+m,[p]))
  probs=np.concatenate((probsneg,probspos,[nu0]))
  u=np.log(h[ys]/h[p]);small=np.abs(ys-p)<=p/2
  L=u*small
  kval=np.arange(math.ceil(p/2),p-1,dtype=int)
  swallow=.5*np.sum(d[kval+1]/(2*kval+1)*np.sum(w[None,:]/(m[None,:]+(kval+1)[:,None]),axis=1)*h[p-kval-1]/h[p]) if len(kval) else 0.
  delta=1-nu0-np.sum(pos)-np.sum(neg)
  res.append({'model':name,'p':p,'row_error':float(abs(np.sum(probs)-1)),'R':float(rate),'big_over_R':float(np.sum(probs[~small])/rate),'truncated_mean_over_R':float(abs(np.sum(probs*L))/rate),'truncated_second_over_R':float(np.sum(probs*L**2)/rate),'delta_over_R':float(delta/rate),'b_over_R':float(swallow/rate)})
pathlib.Path('universal_checks.json').write_text(json.dumps(res,indent=2),encoding='utf8')
print('rows',len(res),'max normalization error',max(a['row_error'] for a in res))
for key in ['big_over_R','truncated_mean_over_R','truncated_second_over_R','delta_over_R']:
 print(key,max(a[key] for a in res))
print('min b/R for p>=8',min(a['b_over_R'] for a in res if a['p']>=8))
