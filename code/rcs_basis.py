"""Harrell restricted cubic spline basis (3 df) on age with fixed published knots.
KNOTS=[41,53,65,80]; basis columns: age (linear), r1, r2.
Self-contained: the risk calculator uses these exact formulas -> reproducible without patsy.
"""
import numpy as np
KNOTS=[41.0,53.0,65.0,80.0]
def _cp(x,t): return np.maximum(x-t,0.0)**3
def rcs(age):
    t1,t2,t3,t4=KNOTS
    x=np.asarray(age,dtype=float)
    d=t4-t3
    r1=_cp(x,t1)-_cp(x,t3)*(t4-t1)/d+_cp(x,t4)*(t3-t1)/d
    r2=_cp(x,t2)-_cp(x,t3)*(t4-t2)/d+_cp(x,t4)*(t3-t2)/d
    return x,r1,r2
def add_basis(df):
    a,r1,r2=rcs(df['age'].values)
    out=df.copy(); out['age_lin']=a; out['age_r1']=r1; out['age_r2']=r2
    return out
if __name__=="__main__":
    print("rcs(70)=",rcs([70.0]))
