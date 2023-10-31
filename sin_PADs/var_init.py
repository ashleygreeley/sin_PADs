import numpy as np
import sin_PADs.fit_funcs

def initiate_variables(L3,instrument,particles,rel):
    B_eq=L3['B_Eq'][:]
    B_calc=L3['B_Calc'][:]

    if ((instrument=='rept')&(particles=='electrons')):
        print('rept electrons')
        L3time=L3['Epoch'][:]
        L3alpha=L3['FEDU_Unbinned_Alpha'][:] #[7000,36]
        L3alpha_binned=L3['FEDU_Alpha'][:]
        L3alphaindex=L3['FEDU_Unbinned_Alpha'][0,:]
        L3fedu=L3['FEDU_Unbinned_0to180']
        L3energy=L3['FEDU_Energy'][:]
        return L3time,L3alpha,L3alpha_binned,L3alphaindex,L3fedu,L3energy,B_eq,B_calc
    if ((instrument=='mageis')&(particles=='electrons')):
        print('mageis electrons')
        L3time=L3['Epoch'][:]
        L3alpha=L3['FEDU_Alpha'][:] #[11]
        L3alpha_binned=L3['FEDU_Alpha'][:] #[11]
        L3alphaindex=L3['FEDU_Alpha'][:]
        L3fedu=L3['FEDU']
        if rel=='rel04':
            L3energy=L3['FEDU_Energy'][:] #1820
        elif rel=='rel03':
            L3energy=L3['FEDU_Energy'][0,:] #1820
        return L3time, L3alpha,L3alpha_binned,L3alphaindex,L3fedu,L3energy,B_eq,B_calc
    if ((instrument=='mageis')&(particles=='protons')):
        print('mageis protons')
        L3time=L3['Epoch'][:]
        L3alpha=L3['FPDU_Alpha'][:]
        L3alpha_binned=L3['FPDU_Alpha'][:]
        L3alphaindex=L3['FPDU_Alpha'][:]
        L3fpdu=L3['FPDU'][:]
        L3energy=L3['FPDU_Energy'][:]
        if rel=='rel04':
            L3energy=L3['FPDU_Energy'][:] #1820
        elif rel=='rel03':
            L3energy=L3['FPDU_Energy'][0,:] #1820
        return L3time,L3alpha,L3alpha_binned,L3alphaindex,L3fpdu,L3energy,B_eq,B_calc
    if ((instrument=='rept')&(particles=='protons')):
        print('rept protons')
        L3time=L3['Epoch'][:]
        L3alpha=L3['FPDU_Unbinned_Alpha'][:] #[7000,36]
        L3alpha_binned=L3['FPDU_Alpha'][:] #[7000,36]
        L3alphaindex=L3['FPDU_Unbinned_Alpha'][0,:]
        L3fpdu=L3['FPDU_Unbinned_0to180']
        L3energy=L3['FPDU_Energy'][0,:]
        return L3time,L3alpha,L3alpha_binned,L3alphaindex,L3fpdu,L3energy,B_eq,B_calc
def create_arrays(time_len,alpha_len,energy_len,legendre):
    flux_fromfit_arr=np.nan*np.ones((time_len,alpha_len,energy_len), dtype=float)
    flux_fromfit_arr_atbinned=np.nan*np.ones((time_len,alpha_len,energy_len), dtype=float)
    rmse_arr=np.nan*np.ones((time_len,energy_len), dtype=float);
    butt_arr=np.nan*np.ones((time_len,energy_len), dtype=float);
    PAIndex_arr=np.nan*np.ones((time_len,energy_len), dtype=float);
    fit_arr=np.nan*np.ones((time_len,energy_len), dtype=float);
    fit_arr_atbinned=np.nan*np.ones((time_len,energy_len), dtype=float);
    maxloarr=np.nan*np.ones((time_len,energy_len), dtype=float);
    if legendre=='yes':
        cn=np.nan*np.ones((time_len,energy_len,5), dtype=float);
        cfull=np.nan*np.ones((time_len,energy_len,5), dtype=float);
        l_rmse_arr=np.nan*np.ones((time_len,energy_len), dtype=float);
        l_fit_arr=np.nan*np.ones((time_len,alpha_len,energy_len), dtype=float);
        legendre_arrs=[cn,cfull,l_rmse_arr,l_fit_arr]
        #return flux_fromfit_arr,flux_fromfit_arr_atbinned,rmse_arr,butt_arr,PAIndex_arr,fit_arr,fit_arr_atbinned,maxloarr,cn,cfull,l_rmse_arr,l_fit_arr
    else:
        legendre_arrs=[]
        #return flux_fromfit_arr,flux_fromfit_arr_atbinned,rmse_arr,butt_arr,PAIndex_arr,fit_arr,fit_arr_atbinned,maxloarr
    return flux_fromfit_arr,flux_fromfit_arr_atbinned,rmse_arr,butt_arr,PAIndex_arr,fit_arr,fit_arr_atbinned,maxloarr,legendre_arrs
def all_fits(tim,en,eq_pitch,L3fedu,L3alpha,L3alpha_binned,legendre):
        xfit1=eq_pitch[tim]
        #7790,11,25
        yfit1=L3fedu[tim,:,en]
        """
        Mirror around 90 degrees
        Important for mageis data
        """

        xfit2=np.abs(xfit1-180)
        yfit2=yfit1

        xfittot=np.append(xfit1,xfit2)
        yfittot=np.append(yfit1,yfit2)
        #print(yfittot)

        fitparams_eq=sin_PADs.fit_funcs.curvetest(xfittot,yfittot);
        maxlo=max(xfittot[xfittot<90]);

        yfit=sin_PADs.fit_funcs.test_func(xfit1*np.pi/180.,fitparams_eq[0],fitparams_eq[1]);
        yfit_atbinned=sin_PADs.fit_funcs.test_func(L3alpha_binned*np.pi/180.,fitparams_eq[0],fitparams_eq[1]);
        rmse=sin_PADs.fit_funcs.RMSEtestnonint(yfit,yfit1);
        butterfly=sin_PADs.fit_funcs.butterflytest(xfit1,yfit1);

        if legendre=='yes':
            if len(xfittot[(yfittot>0)&(xfittot>0)]) > 3:
                cc=np.polynomial.legendre.legfit(np.cos(np.array(xfittot[(yfittot>0)&(xfittot>0)])*np.pi/180.),yfittot[(yfittot>0)&(xfittot>0)],4)
                cfull[tim,en]=cc
                cc[1:]=cc[1:]/cc[0]
            else:
                cc=[np.nan,np.nan,np.nan,np.nan,np.nan]
                cfull[tim,en]=cc

            cn[tim,en]=cc

            cs=np.cos(np.array(xfit1)*np.pi/180.)
            l_fit=cfull[tim,en,0]*legendre(0)(cs)+cfull[tim,en,1]*legendre(1)(cs)+cfull[tim,en,2]*legendre(2)(cs)+cfull[tim,en,3]*legendre(3)(cs)+cfull[tim,en,4]*legendre(4)(cs)
            l_rmse=sin_PADs.fit_funcs.RMSEtestnonint(np.log10(yfit1),l_fit);
            return yfit,yfit_atbinned,rmse,butterfly,fitparams_eq,maxlo,l_fit,l_rmse,cn,cs
        else:
            return yfit,yfit_atbinned,rmse,butterfly,fitparams_eq,maxlo
