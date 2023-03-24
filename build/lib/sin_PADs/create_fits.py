import os
#os.environ["CDF_LIB"] ="/home/jovyan/efs/agreeley/cdf38_1-dist/lib"
from spacepy import pycdf
import numpy as np
import math
from scipy import optimize
import glob
import h5py
import datetime
import matplotlib.pyplot as plt
import datetime
import pathlib
from os.path import exists
from scipy.special import legendre

import sin_PADs
from sin_PADs import config
import sin_PADs.fit_funcs
from sin_PADs.read_ephem import readephemdata
#sudo python setup.py install
#sudo python -m sin_PADs config
def load_data(sc_id, day,instrument,rewrite):
    """

    Parameters
    ----------
    sc_id: str
        The spacecraft id, either "A" or "B". Case insensitive.
    day: str or sin_PADs.create_fits.load_data('B',date,'mageis','yes')datetime.datetime
        A string specifying the full date YYYY-MM-DD
    rewrite: if cdf already exists, rewrite? yes or no

    Returns
    -------

    Example
    -------
    sc_id = 'A'
    day = '2020-01-01'
    """

    ephemDir=sin_PADs.project_data_dir+'/RBSP'
    RBSPDir=sin_PADs.project_data_dir+'/RBSP'

    d1=datetime.datetime.strptime(day,'%Y-%m-%d')
    yr=d1.year
    mo=f'{d1.month:02d}'
    day=f'{d1.day:02d}'


    print(RBSPDir)

    if instrument=='rept':
        print(RBSPDir+'/rept/L2/'+sc_id+'/rbsp'+sc_id.lower()+'_rel03_ect-rept-sci-*2_'+str(yr)+mo+day+'_v*.cdf')
        L2files=glob.glob(RBSPDir+'/rept/L2/'+sc_id+'/rbsp'+sc_id.lower()+'_rel03_ect-rept-sci-*2_'+str(yr)+mo+day+'_v*.cdf')
        L3files=glob.glob(RBSPDir+'/rept/L3/'+sc_id+'/rbsp'+sc_id.lower()+'_rel03_ect-rept-sci-*3_'+str(yr)+mo+day+'_v*.cdf')
    if instrument=='mageis':
        L2files=glob.glob(RBSPDir+'/mageis/L2/'+sc_id+'/rbsp'+sc_id.lower()+'_rel0*_ect-mageis-*2_'+str(yr)+str(mo)+str(day)+'_v*.cdf')
        L3files=glob.glob(RBSPDir+'/mageis/L3/'+sc_id+'/rbsp'+sc_id.lower()+'_rel0*_ect-mageis-*3_'+str(yr)+mo+day+'_v*.cdf')
    L2 = pycdf.CDF(L2files[0])
    L3 = pycdf.CDF(L3files[0])

    print(sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id+'/rbsp'+sc_id.lower()+'_ect_'+instrument+'-PAI'+str(yr)+str(mo)+str(day)+'.cdf')

    file_exists=exists(sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id+'/rbsp'+sc_id.lower()+'_ect_'+instrument+'-PAI'+str(yr)+str(mo)+str(day)+'.cdf')

    if ((file_exists)&(rewrite=='yes')):
        print('overwriting cdf')
        os.remove(sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id+'/rbsp'+sc_id.lower()+'_ect_'+instrument+'-PAI'+str(yr)+str(mo)+str(day)+'.cdf')
        fit_data(L2,L3,sc_id,yr,mo,day,instrument)
    if ((file_exists==False)):
        fit_data(L2,L3,sc_id,yr,mo,day,instrument)
    if ((file_exists)&(rewrite=='no')):
        print('Skipping this one')
        pass


def fit_data(L2,L3,sc_id,yr,mo,day,instrument):

    xxA=[pyind.timetuple().tm_yday+float(pyind.strftime("%H"))/24.+float(pyind.strftime("%M"))/1400.+float(pyind.strftime("%S"))/86400. for pyind in L3['Epoch']]
    ephem=readephemdata(sc_id.lower(),yr,mo,day)
    xx2A=ephem[1]
    orbnum = np.searchsorted(xx2A,xxA, side="left")
    inout=ephem[0]

    if instrument=='rept':
        L3time=L3['Epoch'][:]
        L3alpha=L3['FEDU_Unbinned_Alpha'][:] #[7000,36]
        L3alphaindex=L3['FEDU_Unbinned_Alpha'][0,:]
        L3fedu=L3['FEDU_Unbinned_0to180']
        L3energy=L3['FEDU_Energy'][:]
    if instrument=='mageis':
        L3time=L3['Epoch'][:]
        L3alpha=L3['FEDU_Alpha'][:] #[11]
        L3alphaindex=L3['FEDU_Alpha'][:]
        L3fedu=L3['FEDU']
        L3energy=L3['FEDU_Energy'][:] #1820


    eq_pitch=sin_PADs.fit_funcs.pitchangle_func(L3alpha.astype('float'),L3['B_Eq'][:],L3['B_Calc'][:],instrument)
    eq_pitch_binned=sin_PADs.fit_funcs.pitchangle_func(L3['FEDU_Alpha'][:],L3['B_Eq'][:],L3['B_Calc'][:],instrument,'1d')

    flux_fromfit_arr=np.nan*np.ones((len(L3time),len(L3alphaindex),len(L3energy)), dtype=float)
    flux_fromfit_arr_atbinned=np.nan*np.ones((len(L3time),len(L3['FEDU_Alpha'][:]),len(L3energy)), dtype=float)
    rmse_arr=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);
    l_rmse_arr=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);
    l_fit_arr=np.nan*np.ones((len(L3time),len(L3alphaindex),len(L3energy)), dtype=float);
    butt_arr=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);
    PAIndex_arr=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);
    fit_arr=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);
    fit_arr_atbinned=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);
    cn=np.nan*np.ones((len(L3time),len(L3energy),5), dtype=float);
    cfull=np.nan*np.ones((len(L3time),len(L3energy),5), dtype=float);
    maxloarr=np.nan*np.ones((len(L3time),len(L3energy)), dtype=float);



    for en in range(0,len(L3energy[:])):
        for tim in range(0,len(L3time)):

            #assume 0-90 = 90-180
            xfit1=eq_pitch[tim]
            yfit1=L3fedu[tim,:,en]

            xfit2=np.abs(xfit1-180)
            yfit2=yfit1

            xfittot=np.append(xfit1,xfit2)
            yfittot=np.append(yfit1,yfit2)


            fitparams_eq=sin_PADs.fit_funcs.curvetest(xfittot,yfittot);
            #closest PA to 90
            maxlo=max(xfittot[xfittot<90]);


            yfit=sin_PADs.fit_funcs.test_func(xfit1*np.pi/180.,fitparams_eq[0],fitparams_eq[1]);
            yfit_atbinned=sin_PADs.fit_funcs.test_func(L3['FEDU_Alpha'][:]*np.pi/180.,fitparams_eq[0],fitparams_eq[1]);
            rmse=sin_PADs.fit_funcs.RMSEtestnonint(yfit,yfit1);
            butterfly=sin_PADs.fit_funcs.butterflytest(xfit1,yfit);
            #a*sin**B
            #a
            fit_arr[tim,en]=fitparams_eq[0]
            #b
            PAIndex_arr[tim,en]=fitparams_eq[1]

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


            for alph in range(0,11):
                flux_fromfit_arr_atbinned[tim,alph,en]=yfit_atbinned[alph];
            for alph in range(0,len(yfit)):
                #y values from the fit?
                flux_fromfit_arr[tim,alph,en]=yfit[alph];

                l_fit_arr[tim,alph,en]=l_fit[alph]
            rmse_arr[tim,en]=rmse;
            butt_arr[tim,en]=butterfly;
            l_rmse_arr[tim,en]=l_rmse
            maxloarr[tim,en]=maxlo

    print('Writing to cdf')
    cdf2 = pycdf.CDF(sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id+'/rbsp'+sc_id.lower()+'_ect_'+instrument+'-PAI'+str(yr)+str(mo)+str(day)+'.cdf', '')

    #Depend, time
    cdf2['Epoch']=L3['Epoch'][:]
    cdf2['Epoch'].attrs['SCALETYP']='linear'
    cdf2['Epoch'].attrs['LABLAXIS']='Time'

    #Depend, energy
    #keV for mageis, MeV for rept
    cdf2['Energy']=L3['FEDU_Energy'][:]
    cdf2['Energy'].attrs['FIELDNAM']='Energy'
    cdf2['Energy'].attrs['CATDESC']='Energy bins from FEDU'
    cdf2['Energy'].attrs['FILLVAL']=-1e+31
    cdf2['Energy'].attrs['VALIDMIN']=0.0
    cdf2['Energy'].attrs['VALIXMAX']=100.0
    cdf2['Energy'].attrs['FORMAT']='E12.2'
    cdf2['Energy'].attrs['LABLAXIS']='Energy'
    cdf2['Energy'].attrs['SCALETYP']='log'



    if instrument=='mageis':
        cdf2['EnergyMeV']=L3['FEDU_Energy'][:]/1000.
        cdf2['EnergyMeV'].attrs['UNITS']='MeV'
        cdf2['EnergyMeV'].attrs['VALIDMIN']=0.0
        cdf2['EnergyMeV'].attrs['VALIXMAX']=100.0
        cdf2['EnergyMeV'].attrs['LABLAXIS']='Energy (MeV)'

        cdf2['FESA_CORR']=L2['FESA_CORR'][:]*1000.
        cdf2['FESA_CORR'].attrs['DEPEND_0']= 'Epoch'
        cdf2['FESA_CORR'].attrs['DEPEND_1']= 'EnergyMeV'
        cdf2['FESA_CORR'].attrs['VAR_TYPE']='data'





    #for unbinned fedu 0 to 180
    cdf2['Alpha_equatorial']=eq_pitch[:]
    cdf2['Alpha_equatorial'].attrs['VALIDMIN']=0.0
    cdf2['Alpha_equatorial'].attrs['VALIDMAX']=180.0
    cdf2['Alpha_equatorial'].attrs['LABLAXIS']='Equatorial Pitch Angle'
    cdf2['Alpha_equatorial'].attrs['FIELDNAM']='Eq_Alpha'
    cdf2['Alpha_equatorial'].attrs['UNITS']='degrees'
    cdf2['Alpha_equatorial'].attrs['SCALETYP']='linear'
    cdf2['Alpha_equatorial'].attrs['DEPEND_0']= 'Epoch'
    cdf2['Alpha_equatorial'].attrs['CATDESC']='Pitch Angle mapped to equator'
    cdf2['Alpha_equatorial'].attrs['VAR_NOTES']='FEDU_Unbinned_Alpha to equatorial'
    cdf2['Alpha_equatorial'].attrs['VAR_TYPE']='support_data'

    #for binned fedu 0 to 180
    cdf2['Alpha_equatorial_binned']=eq_pitch_binned[:]
    cdf2['Alpha_equatorial_binned'].attrs['VALIDMIN']=0.0
    cdf2['Alpha_equatorial_binned'].attrs['VALIDMAX']=180.0
    cdf2['Alpha_equatorial_binned'].attrs['LABLAXIS']='Equatorial Pitch Angle'
    cdf2['Alpha_equatorial_binned'].attrs['FIELDNAM']='Eq_Alpha_Binned'
    cdf2['Alpha_equatorial_binned'].attrs['UNITS']='degrees'
    cdf2['Alpha_equatorial_binned'].attrs['SCALETYP']='linear'
    cdf2['Alpha_equatorial_binned'].attrs['DEPEND_0']= 'Epoch'
    cdf2['Alpha_equatorial_binned'].attrs['CATDESC']='Pitch Angle mapped to equator'
    cdf2['Alpha_equatorial_binned'].attrs['VAR_NOTES']='FEDU_Alpha to equatorial'
    cdf2['Alpha_equatorial_binned'].attrs['VAR_TYPE']='support_data'

    cdf2['FEDU_Alpha']=L3['FEDU_Alpha'][:]
    cdf2['FEDU_Alpha'].attrs['VALIDMIN']=0.0
    cdf2['FEDU_Alpha'].attrs['VALIDMAX']=180.0
    cdf2['FEDU_Alpha'].attrs['FIELDNAM']='FEDU_Alpha'
    cdf2['FEDU_Alpha'].attrs['CATDESC']='Pitch Angles for FEDU, binned'
    cdf2['FEDU_Alpha'].attrs['VAR_TYPE']='support_data'


    if instrument=='rept':
        cdf2['FEDU_Unbinned_Alpha']=L3alpha
        cdf2['FEDU_Unbinned_Alpha'].attrs['VALIDMIN']=0.0
        cdf2['FEDU_Unbinned_Alpha'].attrs['VALIDMAX']=180.0
        cdf2['FEDU_Unbinned_Alpha'].attrs['FIELDNAM']='FEDU_Alpha'

    cdf2['flux_fromfit_eq']=flux_fromfit_arr
    cdf2['flux_fromfit_eq'].attrs['units'] = 'cm2sr'
    cdf2['flux_fromfit_eq'].attrs['DEPEND_0']= 'Epoch'
    cdf2['flux_fromfit_eq'].attrs['DEPEND_1']= 'Alpha_equatorial'
    cdf2['flux_fromfit_eq'].attrs['VAR_NOTES']='flux values from the fitting algorithm (to compare data with fit).'
    if instrument=='mageis':
        cdf2['flux_fromfit_eq'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['flux_fromfit_eq'].attrs['DEPEND_2']= 'Energy'

    #L3alpha=L3['FEDU_Unbinned_Alpha'][:] #[7000,36]
    #cdf2['flux_fit']=flux_fromfit_arr
    #cdf2['flux_fit'].attrs['units'] = 'cm2sr'
    #cdf2['flux_fit'].attrs['DEPEND_0']= 'Epoch'
    #cdf2['flux_fit'].attrs['DEPEND_1']= 'FEDU_Alpha'
    #if instrument=='mageis':
    #    cdf2['flux_fit'].attrs['DEPEND_2']= 'EnergyMeV'
    #if instrument=='rept':
    #    cdf2['flux_fit'].attrs['DEPEND_2']= 'Energy'

    cdf2['flux_fromfit_binned']=flux_fromfit_arr_atbinned
    cdf2['flux_fromfit_binned'].attrs['units'] = 'cm2sr'
    cdf2['flux_fromfit_binned'].attrs['DEPEND_0']= 'Epoch'
    cdf2['flux_fromfit_binned'].attrs['DEPEND_1']= 'FEDU_Alpha'
    cdf2['flux_fromfit_binned'].attrs['VAR_NOTES']='flux values from the fitting algorithm (to compare data with fit).'
    if instrument=='mageis':
        cdf2['flux_fromfit_binned'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['flux_fromfit_binned'].attrs['DEPEND_2']= 'Energy'
    #calc_flux_fromfit
    #flux_from_fit
    cdf2['flux_fromfit_legendre']=l_fit_arr
    cdf2['flux_fromfit_legendre'].attrs['units'] = 'cm2sr'
    cdf2['flux_fromfit_legendre'].attrs['DEPEND_0']= 'Epoch'
    cdf2['flux_fromfit_legendre'].attrs['DEPEND_1']= 'FEDU_Unbinned_Alpha'
    cdf2['flux_fromfit_legendre'].attrs['VAR_NOTES']='flux values from the fitting algorithm (to compare data with fit).'
    if instrument=='mageis':
        cdf2['flux_fromfit_legendre'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['flux_fromfit_legendre'].attrs['DEPEND_2']= 'Energy'

    #PAI, time and energy depend
    cdf2['PAIndex']=PAIndex_arr
    cdf2['PAIndex'].attrs['DEPEND_0']= 'Epoch'
    cdf2['PAIndex'].attrs['VALIDMIN']=0
    cdf2['PAIndex'].attrs['VALIDMAX']=10
    cdf2['PAIndex'].attrs['VAR_TYPE']='data'
    cdf2['PAIndex'].attrs['DISPLAY_TYPE']='time_series'
    if instrument=='mageis':
        cdf2['PAIndex'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['PAIndex'].attrs['DEPEND_1']= 'Energy'

    #legendre polynomials
    cdf2['cn']=cn
    cdf2['cn'].attrs['DEPEND_0']= 'Epoch'
    cdf2['cn'].attrs['DISPLAY_TYPE']='time_series'
    cdf2['cn'].attrs['VAR_TYPE']='data'
    if instrument=='mageis':
        cdf2['cn'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['cn'].attrs['DEPEND_1']= 'Energy'

    cdf2['cn2b4']=cn[:,:,2]/cn[:,:,4]
    cdf2['cn2b4'].attrs['DEPEND_0']= 'Epoch'
    cdf2['cn2b4'].attrs['DISPLAY_TYPE']='time_series'
    cdf2['cn2b4'].attrs['VAR_TYPE']='data'
    if instrument=='mageis':
        cdf2['cn2b4'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['cn2b4'].attrs['DEPEND_1']= 'Energy'

    cdf2['cfull']=cfull
    cdf2['cfull'].attrs['DEPEND_0']= 'Epoch'
    cdf2['cfull'].attrs['DISPLAY_TYPE']='time_series'
    cdf2['cfull'].attrs['VAR_TYPE']='data'
    if instrument=='mageis':
        cdf2['cfull'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['cfull'].attrs['DEPEND_1']= 'Energy'


    cdf2['butterfly']=butt_arr
    cdf2['butterfly'].attrs['DEPEND_0']= 'Epoch'
    cdf2['butterfly'].attrs['FIELDNAM']='butterfly'
    cdf2['butterfly'].attrs['VAR_TYPE']='support_data'
    cdf2['butterfly'].attrs['DISPLAY_TYPE']='time_series'
    cdf2['butterfly'].attrs['CATDESC']='binary butterfly from sin fit'
    cdf2['butterfly'].attrs['VAR_NOTES']='0 for non butterfly, 1 for butterfly PAD'
    if instrument=='mageis':
        cdf2['butterfly'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['butterfly'].attrs['DEPEND_1']= 'Energy'



    cdf2['J0']=fit_arr
    cdf2['J0'].attrs['DEPEND_0']= 'Epoch'
    if instrument=='mageis':
        cdf2['J0'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['J0'].attrs['DEPEND_1']= 'Energy'
    cdf2['J0'].attrs['VAR_TYPE']='support_data'
    cdf2['J0'].attrs['CATDESC']='parameter 0 (a) from fit a*sin(pa)**b'
    cdf2['J0'].attrs['DISPLAY_TYPE']='time_series'


    cdf2['RMSE_legendre']=l_rmse_arr
    cdf2['RMSE_legendre'].attrs['DEPEND_0']= 'Epoch'
    cdf2['RMSE_legendre'].attrs['VALIDMIN']=0
    cdf2['RMSE_legendre'].attrs['VALIDMAX']=1
    cdf2['RMSE_legendre'].attrs['FIELDNAM']='RMSE'
    cdf2['RMSE_legendre'].attrs['VAR_TYPE']='support_data'
    if instrument=='mageis':
        cdf2['RMSE_legendre'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['RMSE_legendre'].attrs['DEPEND_1']= 'Energy'
    cdf2['RMSE_legendre'].attrs['CATDESC']='RMSE for legendre fit'


    cdf2['RMSE']=rmse_arr
    cdf2['RMSE'].attrs['DEPEND_0']= 'Epoch'
    cdf2['RMSE'].attrs['VALIDMIN']=0
    cdf2['RMSE'].attrs['VALIDMAX']=1
    cdf2['RMSE'].attrs['FIELDNAM']='RMSE'
    cdf2['RMSE'].attrs['VAR_TYPE']='support_data'
    if instrument=='mageis':
        cdf2['RMSE'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['RMSE'].attrs['DEPEND_1']= 'Energy'

    if instrument=='mageis':
        cdf2['FESA']=L2['FESA'][:]*1000.
        cdf2['FESA'].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['FESA']=L2['FESA'][:]
        cdf2['FESA'].attrs['DEPEND_1']= 'Energy'

    cdf2['FESA'].attrs['DEPEND_0']= 'Epoch'
    cdf2['FESA'].attrs['SCALETYP'] = 'log'
    cdf2['FESA'].attrs['VALIDMIN'] = 0.01
    cdf2['FESA'].attrs['VALIDMAX']=1e+31
    cdf2['FESA'].attrs['CATDESC']='Spin-averaged differential electron flux (lowest electron bin is not valid inside of L=2.8) and integral flux for the last 2 bins.'
    #cdf2['FESA'].attrs['VAR_NOTES']='MeV for mageis and rept, not corr for MagEIS'
    cdf2['FESA'].attrs['DISPLAY_TYPE']='nnSpectrogram'
    cdf2['FESA'].attrs['FIELDNAM']='FESA'
    cdf2['FESA'].attrs['FILLVAL']=-1e+31
    cdf2['FESA'].attrs['FORMAT']='E10.3'
    cdf2['FESA'].attrs['VAR_TYPE']='data'

    if instrument=='mageis':
        cdf2['FEDU']=L3['FEDU'][:]*1000.
        cdf2['FEDU'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['FEDU']=L3['FEDU'][:]
        cdf2['FEDU'].attrs['DEPEND_2']= 'Energy'
    cdf2['FEDU'].attrs['DEPEND_0']= 'Epoch'
    cdf2['FEDU'].attrs['DEPEND_1']= 'FEDU_Alpha' #FEDU_Alpha is length 17
    cdf2['FEDU'].attrs['DISPLAY_TYPE']='nnSpectrogram'
    cdf2['FEDU'].attrs['FIELDNAM']='FEDU'
    cdf2['FEDU'].attrs['FILLVAL']=-1e+31
    cdf2['FEDU'].attrs['FORMAT']='E12.2'
    cdf2['FEDU'].attrs['SCALTYP']='log'
    cdf2['FEDU'].attrs['VAR_TYPE']='data'
    cdf2['FEDU'].attrs['VALIDMIN']=0
    cdf2['FEDU'].attrs['VALIDMAX']=1e+31

    #binned FEDU with equatorial pitch angles
    if instrument=='mageis':
        cdf2['FEDU_Eq']=L3['FEDU'][:]*1000.
        cdf2['FEDU_Eq'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['FEDU_Eq']=L3['FEDU'][:]
        cdf2['FEDU_Eq'].attrs['DEPEND_2']= 'Energy'
    cdf2['FEDU_Eq'].attrs['DEPEND_0']= 'Epoch'
    cdf2['FEDU_Eq'].attrs['DEPEND_1']= 'Alpha_equatorial_binned'
    cdf2['FEDU_Eq'].attrs['DISPLAY_TYPE']='nnSpectrogram'
    cdf2['FEDU_Eq'].attrs['FIELDNAM']='FEDU_Eq'
    cdf2['FEDU_Eq'].attrs['FILLVAL']=-1e+31
    cdf2['FEDU_Eq'].attrs['FORMAT']='E12.2'
    cdf2['FEDU_Eq'].attrs['SCALTYP']='log'
    cdf2['FEDU_Eq'].attrs['VAR_TYPE']='data'
    cdf2['FEDU_Eq'].attrs['VALIDMIN']=0
    cdf2['FEDU_Eq'].attrs['VALIDMAX']=1e+31

    cdf2['MLAT']=L3['MLAT'][:]
    cdf2['MLAT'].attrs['DEPEND_0']= 'Epoch'
    cdf2['MLAT'].attrs['FIELDNAM']='MLAT'
    cdf2['MLAT'].attrs['FILLVAL']=-1e+31
    cdf2['MLAT'].attrs['FORMAT']='E10.3'
    cdf2['MLAT'].attrs['LABLAXIS']='MLAT'
    cdf2['MLAT'].attrs['SCALETYP']='linear'
    cdf2['MLAT'].attrs['UNITS']='degrees'
    cdf2['MLAT'].attrs['VALIDMAX']=90.0
    cdf2['MLAT'].attrs['VALIDMIN']=-90.0
    cdf2['MLAT'].attrs['VAR_TYPE']='support_data'
    cdf2['MLAT'].attrs['CATDESC']='MLAT - Eccentric Dipole Magnetic Latitude. Computed using OP77Q external field and IGRF internal field.'
    cdf2['MLAT'].attrs['DISPLAY_TYPE']='series'

    cdf2['MLT']=L3['MLT'][:]
    cdf2['MLT'].attrs['DEPEND_0']= 'Epoch'
    cdf2['MLT'].attrs['DISPLAY_TYPE']='time_series'
    cdf2['MLT'].attrs['LABLAXIS']='MLT'
    cdf2['MLT'].attrs['SI_conversion']='2.778e-4>s'
    cdf2['MLT'].attrs['UNITS']='h'
    cdf2['MLT'].attrs['VAR_TYPE']='support_data'
    #cdf2['butterfly']=butt_arr

    cdf2['Lstar']=L3['L_star'][:]
    cdf2['Lstar'].attrs['DEPEND_0']= 'Epoch'
    cdf2['Lstar'].attrs['CATDESC']='Lstar'
    cdf2['Lstar'].attrs['LABLAXIS']='L*'
    cdf2['Lstar'].attrs['FIELDNAM']='Lstar'
    cdf2['Lstar'].attrs['SCALETYP']='linear'
    cdf2['Lstar'].attrs['VALIDMIN']=0
    cdf2['Lstar'].attrs['VALIDMAX']=30
    cdf2['Lstar'].attrs['VAR_TYPE']='support_data'

    cdf2['L']=L3['L'][:]
    cdf2['L'].attrs['DEPEND_0']= 'Epoch'
    cdf2['L'].attrs['CATDESC']='L'
    cdf2['L'].attrs['LABLAXIS']='L'
    cdf2['L'].attrs['FIELDNAM']='L'
    cdf2['L'].attrs['SCALETYP']='linear'
    cdf2['L'].attrs['VALIDMIN']=0
    cdf2['L'].attrs['VALIDMAX']=30
    cdf2['L'].attrs['VAR_TYPE']='support_data'

    cdf2['B_Calc']=L2['B_Calc'][:]
    cdf2['B_Calc'].attrs['UNITS']='nT'
    cdf2['B_Calc'].attrs['VALIDMIN']=0.0
    cdf2['B_Calc'].attrs['VALIDMAX']=1e+31
    cdf2['B_Calc'].attrs['SCALETYP']='linear'
    cdf2['B_Calc'].attrs['SI_conversion']='1.0e9>T'
    cdf2['B_Calc'].attrs['LABLAXIS']='B'
    cdf2['B_Calc'].attrs['VAR_TYPE']='support_data'

    cdf2['B_Eq']=L2['B_Eq'][:]
    cdf2['B_Eq'].attrs['UNITS']='nT'
    cdf2['B_Eq'].attrs['VALIDMAX']=1e+31
    cdf2['B_Eq'].attrs['SCALETYP']='linear'
    cdf2['B_Eq'].attrs['SI_conversion']='1.0e9>T'
    cdf2['B_Eq'].attrs['LABLAXIS']='B'
    cdf2['B_Eq'].attrs['VAR_TYPE']='support_data'

    print(type(orbnum[0]))
    print(np.shape(orbnum))
    print(L3['Epoch'][0])
    oo=[np.uint8(ii) for ii in orbnum]
    print(np.shape(L3['Epoch'][:]))
    cdf2['orbit_number']=oo
    cdf2['orbit_number'].attrs['DEPEND_0']= 'Epoch'
    cdf2['orbit_number'].attrs['VAR_TYPE']='support_data'


    cdf2['closest_PA_90deg']=maxloarr
    cdf2['closest_PA_90deg'].attrs['DEPEND_0']= 'Epoch'
    cdf2['closest_PA_90deg'].attrs['DEPEND_1']= 'Energy'
    cdf2['closest_PA_90deg'].attrs['VAR_TYPE']='support_data'
    cdf2['closest_PA_90deg'].attrs['CATDESC']='closest PA to 90 degrees'


    if instrument=='rept':
        cdf2['FEDU_Unbinned_0to180']=L3['FEDU_Unbinned_0to180'][:]
        cdf2['FESA'].attrs['UNITS'] = 'cm!e-2!ns!e-1!nsr!e-1!nMeV!e-1!n'
        cdf2['Energy'].attrs['UNITS']='MeV'
        cdf2['FESA'].attrs['VAR_NOTES']='Dimension 1 corresponds to 12 electron energy bins. The units for the last 2 bins are 1/(cm^2-s-sr).'


    cdf2.attrs['Author'] = 'Ashley Greeley'
    cdf2.attrs['CreateDate'] = datetime.datetime.now()
    cdf2.close()


if __name__ == '__main__':
    sc_id = 'A'
    day = '2019-01-26'
    path = load_state(sc_id, day)
