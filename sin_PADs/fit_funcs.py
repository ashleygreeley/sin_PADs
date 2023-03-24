import numpy as np
import math
from scipy import optimize


def pitchangle_func(y,beq,bcalc,instrument,whichd=None):
    """
    Parameters
    ----------
    y : list
        DESCRIPTION.
    beq : TYPE
        DESCRIPTION.
    bcalc : TYPE
        DESCRIPTION.

    Returns
    -------
    pitch_eqtest : TYPE
        DESCRIPTION.
        one beq and bcalc for each time stamp
        y is a list for each timestamp
    """
    if ((instrument=='rept')&(whichd==None)):
        leny=len(y[1])
    if ((instrument=='mageis')&(whichd==None)):
        leny=len(y)
    if whichd=='1d':
        leny=len(y)
    beqsquare=[beq for i in range(0,leny)]
    beqsquare=np.array(beqsquare).T
    bcalcsquare=[bcalc for i in range(0,leny)]
    bcalcsquare=np.array(bcalcsquare).T
    sinsq= np.sin(y.astype('float')*math.pi/180.)**2
    bbyb=beqsquare.astype('float')/bcalcsquare.astype('float')
    sqrtstuff=np.sqrt(bbyb*sinsq)
    pitch_eqtest=np.arcsin(sqrtstuff)*180./math.pi
    if ((instrument=='rept')&(whichd==None)):
        pitch_eqtest[y>90]=180.-pitch_eqtest[y>90]
    if ((instrument=='mageis')or(whichd=='1d')):
        for ee in range(0,len(bcalc)):
            pitch_eqtest[ee,y>90]=180.-pitch_eqtest[ee,y>90]
    return pitch_eqtest


def pitchangle_func2(self,y,beq,bcalc):
    beqsquare=np.column_stack((beq,beq,beq,beq,beq,beq,beq,beq,beq,beq,beq))
    bcalcsquare=np.column_stack((bcalc,bcalc,bcalc,bcalc,bcalc,bcalc,bcalc,bcalc,bcalc,bcalc,bcalc))


    sinsq= np.sin(y.astype('float')*math.pi/180.)**2
    bbyb=beqsquare.astype('float')/bcalcsquare.astype('float')
    sqrtstuff=np.sqrt(bbyb*sinsq)
    pitch_eqtest=np.arcsin(sqrtstuff)*180./math.pi
    for ee in range(0,len(bcalc)):
        pitch_eqtest[ee,y>90]=180.-pitch_eqtest[ee,y>90]
    return pitch_eqtest

def curvetest(xx,yy):
    xx=np.array(xx)
    yy=np.array(yy)
    xx=xx[~np.isnan(yy)]
    yy=yy[~np.isnan(yy)]
    xx=xx[np.isfinite(yy)]
    yy=yy[np.isfinite(yy)]
    pa=np.where((xx>20)&(xx<160)&(yy>0))
    xx=xx[pa]
    yy=yy[pa]
    ff=np.where(yy>0)
    if (len(yy)>3)&(len(ff[0])>3):
        try:
            params, paramsco=optimize.curve_fit(test_func, xx* np.pi / 180. , yy,p0=[8000, 3], maxfev=15000)
            paramarray=[params[0],params[1]]
        except RuntimeError:
            paramarray=[np.nan,np.nan]
    else:
        paramarray=[np.nan,np.nan]
    return    paramarray

def test_func(x, a, b):
    return a * np.abs(np.sin(x))**b

def butterflytest(xx,yy):
    """
    Author: Ashley Greeley
    Purpose: Determine if a pad is butterfly shaped by comparing the mean
    flux between 85 and 95 degrees to the max of the means in a range of
    angles (75-105, 80-100, etc) time 0.95
    Input: Alpha (xx), flux (yy)
    Output: 'butterfly' or 'norm' string
    """
    jj=[]
    xx=np.array(xx)
    yy=np.array(yy)
    for hh in range(6,45):
        gg=np.where((xx > 90-hh) & (xx < 90+hh))
        jj.append(np.mean(yy[gg]))
    gg=np.where((xx > 75) & (xx < 105))
    middle=np.nanmean(yy[gg])
    edges=np.nanmax(jj)*.95
    if middle < edges:
        return 1
    else:
        return 0
def RMSEtest(yydata,yyfit):
    yyfit2=np.array(yyfit)
    yydata2=np.array(yydata)
    yydata3=yydata2[(~np.isnan(yydata2))]
    yyfit3=yyfit2[(~np.isnan(yyfit2))]
    ff=np.where(yyfit3>0)

    if (len(yyfit3>3)&(len(ff[0])>3)):
        rmse=np.sqrt(np.nansum((yydata3-yyfit3)**2)/len(yyfit3))/(max(yyfit3)-min(yyfit3))
    else:
        rmse=np.nan
    return float(rmse)

def RMSEtestnonint(yydata,yyfit):
    yyfit2=np.array(yyfit)
    yydata2=np.array(yydata)

    yydata3=yydata2[(~np.isnan(yydata2))&(~np.isnan(yyfit2))]
    yyfit3=yyfit2[(~np.isnan(yyfit2))&(~np.isnan(yydata2))]
    if len(yyfit3)>3:
        rmse=np.sqrt(np.nansum((yydata3-yyfit3)**2)/len(yyfit3))/(max(yyfit3)-min(yyfit3))
    else:
        rmse=np.nan
    return float(rmse)
