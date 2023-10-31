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
import os.path
from scipy.special import legendre
import pathlib

import sin_PADs
from sin_PADs import config

from typing import Union
import pathlib
import zipfile

import sin_PADs.fit_funcs
from sin_PADs.read_ephem import readephemdata
from sin_PADs.var_init import initiate_variables
from sin_PADs.var_init import create_arrays
from sin_PADs.var_init import all_fits
#sudo python setup.py install
#sudo python -m sin_PADs config
def load_data(sc_id, day,instrument,rewrite,particle,**kwargs):
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
    release= rel03 or rel04. If none, rel04 assumed
    """

    rel = kwargs.get('rel', None)
    if rel=='rel04':
        rel='rel04'
    elif rel=='rel03':
        rel='rel03'
    else:
        print('release number not in proper form, rel04 assumed for mageis and rel03 for rept')
        if instrument=='mageis':
            rel='rel04'
        if instrument=='rept':
            rel='rel03'

    legendre = kwargs.get('legendre', None)
    if legendre=='yes':
        legendre='yes'
    else:
        legenre='no'

    ephemDir=sin_PADs.project_data_dir+'/RBSP'
    RBSPDir=sin_PADs.project_data_dir+'/RBSP'

    d1=datetime.datetime.strptime(day,'%Y-%m-%d')
    yr=d1.year
    mo=f'{d1.month:02d}'
    day=f'{d1.day:02d}'


    #print(RBSPDir)
    file_matchL2='rbsp'+sc_id.lower()+'_'+rel+'_ect-'+instrument+'*2_'+str(yr)+mo+day+'_v*.cdf'
    file_matchL3='rbsp'+sc_id.lower()+'_'+rel+'_ect-'+instrument+'*3_'+str(yr)+mo+day+'_v*.cdf'
    local_filesL2=glob.glob(RBSPDir+'/'+instrument+'/L2/'+sc_id+'/'+file_matchL2)
    local_filesL3=glob.glob(RBSPDir+'/'+instrument+'/L3/'+sc_id+'/'+file_matchL3)




    if len(local_filesL2) in [1, 2]:
        print('at least on L2 exists in local directory')
        L2 = pycdf.CDF(local_filesL2[0])
    elif len(local_filesL2) == 0:
        print('no local L2, downloading')
        print(f'https://spdf.gsfc.nasa.gov/pub/data/rbsp/rbsp{sc_id.lower()}/l2/ect/{instrument}/sectors/{rel}/{yr}/')
        downloader = sin_PADs.Downloader(
                f'https://spdf.gsfc.nasa.gov/pub/data/rbsp/rbsp{sc_id.lower()}/l2/ect/{instrument}/sectors/{rel}/{yr}/',
                download_dir=RBSPDir+'/'+instrument+'/L2/'+sc_id
                )

        matched_downloaders = downloader.ls(match=file_matchL2)
        file_path = matched_downloaders[0].download(stream='True')
        L2=pycdf.CDF(str(file_path))

    if len(local_filesL3) in [1, 2]:
        print('at least one L3 exists in local directory')
        L3 = pycdf.CDF(local_filesL3[0])
    elif len(local_filesL3) == 0:
        print('no local L3, downloading')
        print(f'https://spdf.gsfc.nasa.gov/pub/data/rbsp/rbsp{sc_id.lower()}/l2/ect/{instrument}/sectors/{rel}/{yr}/')
        downloader = sin_PADs.Downloader(
                f'https://spdf.gsfc.nasa.gov/pub/data/rbsp/rbsp{sc_id.lower()}/l3/ect/{instrument}/sectors/{rel}/{yr}/',
                download_dir=RBSPDir+'/'+instrument+'/L3/'+sc_id
                )

        matched_downloaders = downloader.ls(match=file_matchL3)
        file_path = matched_downloaders[0].download(stream='True')
        L3=pycdf.CDF(str(file_path))


    this_directory=sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id+'/'+particle+'/'
    this_file='rbsp'+sc_id.lower()+'_'+rel+'_ect_'+instrument+'-PAI'+str(yr)+str(mo)+str(day)+'_'+particle+'_v2.cdf'
    this_filename=this_directory+this_file

    print(this_filename)

    file_exists=exists(this_filename)

    if ((file_exists)&(rewrite=='yes')):
        print('overwriting cdf')
        os.remove(this_filename)
        fit_data(L2,L3,sc_id,yr,mo,day,instrument,rel,particle,legendre)
    if ((file_exists==False)):
        print('creating new cdf')
        fit_data(L2,L3,sc_id,yr,mo,day,instrument,rel,particle,legendre)
    if ((file_exists)&(rewrite=='no')):
        print('rewrite=no - Skipping this one')
        pass


def fit_data(L2,L3,sc_id,yr,mo,day,instrument,rel,particle,legendre):

    xxA=[pyind.timetuple().tm_yday+float(pyind.strftime("%H"))/24.+float(pyind.strftime("%M"))/1400.+float(pyind.strftime("%S"))/86400. for pyind in L3['Epoch']]
    ephem=readephemdata(sc_id.lower(),yr,mo,day)
    xx2A=ephem[1]
    orbnum = np.searchsorted(xx2A,xxA, side="left")
    inout=ephem[0]


    if particle=='electrons':
        L3binnedA=L3['FEDU_Alpha'][:]
    elif particle=='protons':
        L3binnedA=L3['FPDU_Alpha'][:]
    L3time,L3alpha,L3alpha_binned,L3alphaindex,L3du,L3energy,B_eq,B_calc=initiate_variables(L3,instrument,particle,rel)


    eq_pitch=sin_PADs.fit_funcs.pitchangle_func(L3alpha.astype('float'),B_eq,B_calc,instrument)
    eq_pitch_binned=sin_PADs.fit_funcs.pitchangle_func(L3binnedA,B_eq,B_calc,instrument,'1d')


    #initialize arrays

    #electron arrays
    #if legendre=='yes':
        #flux_fromfit_arr,flux_fromfit_arr_atbinned,rmse_arr,butt_arr,PAIndex_arr,fit_arr,fit_arr_atbinned,maxloarr,cn,cfull,l_rmse_arr,l_fit_arr=create_arrays(len(L3time),len(L3alphaindex),len(L3energy),legendre)
    #else:

    flux_fromfit_arr,flux_fromfit_arr_atbinned,rmse_arr,butt_arr,PAIndex_arr,fit_arr,fit_arr_atbinned,maxloarr,legendre_arrs=create_arrays(len(L3time),len(L3alphaindex),len(L3energy),legendre)



    for en in range(0,len(L3energy[:])):
        for tim in range(0,len(L3time)):

            #assume 0-90 = 90-180
            #very important for mageis


            yfit,yfit_atbinned,rmse,butterfly,fitparams_eq,maxlo=all_fits(tim,en,eq_pitch,L3du,L3alpha,L3alpha_binned,legendre)


            """
            functional form is y=a*sin(x)**b
            a=fitparams_eq[0]
            b=fitparams_eq[0]
            """
            fit_arr[tim,en]=fitparams_eq[0]
            PAIndex_arr[tim,en]=fitparams_eq[1]


            for alph in range(0,len(L3alpha_binned)):
                flux_fromfit_arr_atbinned[tim,alph,en]=yfit_atbinned[alph];
            for alph in range(0,len(yfit)):
                flux_fromfit_arr[tim,alph,en]=yfit[alph];
                #l_fit_arr[tim,alph,en]=l_fit[alph]
            rmse_arr[tim,en]=rmse;
            butt_arr[tim,en]=butterfly;
            #l_rmse_arr[tim,en]=l_rmse
            maxloarr[tim,en]=maxlo


    print('Writing to cdf')
    filpath=sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id
    if not os.path.isdir(filpath):
        os.makedir(sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id)


    cdf2 = pycdf.CDF(sin_PADs.project_data_dir+'/RBSP/'+instrument+'/L4PAI/'+sc_id+'/'+particle+'/rbsp'+sc_id.lower()+'_'+rel+'_ect_'+instrument+'-PAI'+str(yr)+str(mo)+str(day)+'_'+particle+'_v2.cdf', '')

    #Depend, time
    cdf2['Epoch']=L3['Epoch'][:]
    cdf2['Epoch'].attrs['SCALETYP']='linear'
    cdf2['Epoch'].attrs['LABLAXIS']='Time'


    if particle=='electrons':
        var_name_sa='FESA'
        var_name_du='FEDU'
        var_du_eq='FEDU_Eq'
        du_en=L3['FEDU_Energy'][:]
        du_unbinned='FEDU_Unbinned_0to180'
        var_du_alpha='FEDU_Alpha'
        du_alpha=L3['FEDU_Alpha'][:]
        if instrument=='mageis':
            sa_corr=L2['FESA_CORR'][:]
        if instrument=='rept':
            var_du_unbinned='FEDU_Unbinned_Alpha'
    elif particle=='protons':
        var_name_sa='FPSA'
        var_name_du='FPDU'
        var_du_alpha='FPDU_Alpha'
        du_alpha=L3['FPDU_Alpha'][:]
        du_en=L3['FPDU_Energy'][0,:]
        var_du_eq='FPDU_Eq'
        du_unbinned='FPDU_Unbinned_0to180'
        if instrument=='rept':
            var_du_unbinned='FPDU_Unbinned_Alpha'
    #print('{} {} {}'.format(var_sa,var_du,var_du_alpha))

    #Depend, energy
    #keV for mageis, MeV for rept
    cdf2['Energy']=du_en
    cdf2['Energy'].attrs['FIELDNAM']='Energy'
    cdf2['Energy'].attrs['CATDESC']='Energy bins from '+var_name_du
    cdf2['Energy'].attrs['FILLVAL']=-1e+31
    cdf2['Energy'].attrs['VALIDMIN']=0.0
    cdf2['Energy'].attrs['VALIXMAX']=100.0
    cdf2['Energy'].attrs['FORMAT']='E12.2'
    cdf2['Energy'].attrs['LABLAXIS']='Energy'
    cdf2['Energy'].attrs['SCALETYP']='log'
    if ((instrument=='mageis') & (rel=='rel03')):
        cdf2['Energy'].attrs['DEPEND_0']='Epoch'



    if instrument=='mageis':
        cdf2['EnergyMeV']=du_en/1000.
        cdf2['EnergyMeV'].attrs['UNITS']='MeV'
        cdf2['EnergyMeV'].attrs['VALIDMIN']=0.0
        cdf2['EnergyMeV'].attrs['VALIXMAX']=100.0
        cdf2['EnergyMeV'].attrs['LABLAXIS']='Energy (MeV)'
        if rel=='rel03':
            cdf2['EnergyMeV'].attrs['DEPEND_0']='Epoch'

        if particle=='electrons':
            cdf2['FESA_CORR']=sa_corr*1000.
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
    cdf2['Alpha_equatorial'].attrs['VAR_NOTES']=var_name_du+'_Unbinned_Alpha to equatorial'
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
    cdf2['Alpha_equatorial_binned'].attrs['VAR_NOTES']=var_name_du+'_Alpha to equatorial'
    cdf2['Alpha_equatorial_binned'].attrs['VAR_TYPE']='support_data'





    cdf2[var_du_alpha]=du_alpha
    cdf2[var_du_alpha].attrs['VALIDMIN']=0.0
    cdf2[var_du_alpha].attrs['VALIDMAX']=180.0
    cdf2[var_du_alpha].attrs['FIELDNAM']=var_du_alpha
    cdf2[var_du_alpha].attrs['CATDESC']='Pitch Angles for '+var_name_du+', binned'
    cdf2[var_du_alpha].attrs['VAR_TYPE']='support_data'



    if instrument=='rept':
        cdf2[var_du_unbinned]=L3alpha
        cdf2[var_du_unbinned].attrs['VALIDMIN']=0.0
        cdf2[var_du_unbinned].attrs['VALIDMAX']=180.0
        cdf2[var_du_unbinned].attrs['FIELDNAM']=var_du_alpha

    cdf2['flux_fromfit_eq']=flux_fromfit_arr
    cdf2['flux_fromfit_eq'].attrs['units'] = 'cm2sr'
    cdf2['flux_fromfit_eq'].attrs['DEPEND_0']= 'Epoch'
    cdf2['flux_fromfit_eq'].attrs['DEPEND_1']= 'Alpha_equatorial'
    cdf2['flux_fromfit_eq'].attrs['VAR_NOTES']='flux values from the fitting algorithm (to compare data with fit).'
    if instrument=='mageis':
        cdf2['flux_fromfit_eq'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['flux_fromfit_eq'].attrs['DEPEND_2']= 'Energy'




    cdf2['flux_fromfit_binned']=flux_fromfit_arr_atbinned
    cdf2['flux_fromfit_binned'].attrs['units'] = 'cm2sr'
    cdf2['flux_fromfit_binned'].attrs['DEPEND_0']= 'Epoch'
    cdf2['flux_fromfit_binned'].attrs['DEPEND_1']= var_du_alpha
    cdf2['flux_fromfit_binned'].attrs['VAR_NOTES']='flux values from the fitting algorithm (to compare data with fit).'
    if instrument=='mageis':
        cdf2['flux_fromfit_binned'].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2['flux_fromfit_binned'].attrs['DEPEND_2']= 'Energy'
    #calc_flux_fromfit
    #flux_from_fit




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
    if legendre=='yes':
        cdf2['flux_fromfit_legendre']=l_fit_arr
        cdf2['flux_fromfit_legendre'].attrs['units'] = 'cm2sr'
        cdf2['flux_fromfit_legendre'].attrs['DEPEND_0']= 'Epoch'
        cdf2['flux_fromfit_legendre'].attrs['DEPEND_1']= 'FEDU_Unbinned_Alpha'
        cdf2['flux_fromfit_legendre'].attrs['VAR_NOTES']='flux values from the fitting algorithm (to compare data with fit).'
        if instrument=='mageis':
            cdf2['flux_fromfit_legendre'].attrs['DEPEND_2']= 'EnergyMeV'
        if instrument=='rept':
            cdf2['flux_fromfit_legendre'].attrs['DEPEND_2']= 'Energy'
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








    if instrument=='mageis':
        cdf2[var_name_sa]=L2[var_name_sa][:]*1000.
        cdf2[var_name_sa].attrs['DEPEND_1']= 'EnergyMeV'
    if instrument=='rept':
        cdf2[var_name_sa]=L2[var_name_sa][:]
        cdf2[var_name_sa].attrs['DEPEND_1']= 'Energy'





    cdf2[var_name_sa].attrs['DEPEND_0']= 'Epoch'
    cdf2[var_name_sa].attrs['SCALETYP'] = 'log'
    cdf2[var_name_sa].attrs['VALIDMIN'] = 0.01
    cdf2[var_name_sa].attrs['VALIDMAX']=1e+31
    cdf2[var_name_sa].attrs['CATDESC']='Spin-averaged differential '+particle+' flux.'
    cdf2[var_name_sa].attrs['DISPLAY_TYPE']='nnSpectrogram'
    cdf2[var_name_sa].attrs['FIELDNAM']=var_name_sa
    cdf2[var_name_sa].attrs['FILLVAL']=-1e+31
    cdf2[var_name_sa].attrs['FORMAT']='E10.3'
    cdf2[var_name_sa].attrs['VAR_TYPE']='data'



    if instrument=='mageis':
        cdf2[var_name_du]=L3[var_name_du][:]*1000.
        cdf2[var_name_du].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2[var_name_du]=L3[var_name_du][:]
        cdf2[var_name_du].attrs['DEPEND_2']= 'Energy'
    cdf2[var_name_du].attrs['DEPEND_0']= 'Epoch'
    cdf2[var_name_du].attrs['DEPEND_1']= var_du_alpha #
    cdf2[var_name_du].attrs['DISPLAY_TYPE']='nnSpectrogram'
    cdf2[var_name_du].attrs['FIELDNAM']=var_name_du
    cdf2[var_name_du].attrs['FILLVAL']=-1e+31
    cdf2[var_name_du].attrs['FORMAT']='E12.2'
    cdf2[var_name_du].attrs['SCALTYP']='log'
    cdf2[var_name_du].attrs['VAR_TYPE']='data'
    cdf2[var_name_du].attrs['VALIDMIN']=0
    cdf2[var_name_du].attrs['VALIDMAX']=1e+31



    #binned FEDU with equatorial pitch angles
    if instrument=='mageis':
        cdf2[var_du_eq]=L3[var_name_du][:]*1000.
        cdf2[var_du_eq].attrs['DEPEND_2']= 'EnergyMeV'
    if instrument=='rept':
        cdf2[var_du_eq]=L3['FEDU'][:]
        cdf2[var_du_eq].attrs['DEPEND_2']= 'Energy'
    cdf2[var_du_eq].attrs['DEPEND_0']= 'Epoch'
    cdf2[var_du_eq].attrs['DEPEND_1']= 'Alpha_equatorial_binned'
    cdf2[var_du_eq].attrs['DISPLAY_TYPE']='nnSpectrogram'
    cdf2[var_du_eq].attrs['FIELDNAM']='FEDU_Eq'
    cdf2[var_du_eq].attrs['FILLVAL']=-1e+31
    cdf2[var_du_eq].attrs['FORMAT']='E12.2'
    cdf2[var_du_eq].attrs['SCALTYP']='log'
    cdf2[var_du_eq].attrs['VAR_TYPE']='data'
    cdf2[var_du_eq].attrs['VALIDMIN']=0
    cdf2[var_du_eq].attrs['VALIDMAX']=1e+31

    ###### GENERAL
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

    oo=[np.uint8(ii) for ii in orbnum]
    cdf2['orbit_number']=oo
    cdf2['orbit_number'].attrs['DEPEND_0']= 'Epoch'
    cdf2['orbit_number'].attrs['VAR_TYPE']='support_data'


    cdf2['closest_PA_90deg']=maxloarr
    cdf2['closest_PA_90deg'].attrs['DEPEND_0']= 'Epoch'
    cdf2['closest_PA_90deg'].attrs['DEPEND_1']= 'Energy'
    cdf2['closest_PA_90deg'].attrs['VAR_TYPE']='support_data'
    cdf2['closest_PA_90deg'].attrs['CATDESC']='closest PA to 90 degrees'
    ###########


    if instrument=='rept':
        cdf2[du_unbinned]=L3[du_unbinned][:]
        cdf2[var_name_sa].attrs['UNITS'] = 'cm!e-2!ns!e-1!nsr!e-1!nMeV!e-1!n'
        cdf2['Energy'].attrs['UNITS']='MeV'
        cdf2[var_name_sa].attrs['VAR_NOTES']='Dimension 1 corresponds to 12 electron energy bins. The units for the last 2 bins are 1/(cm^2-s-sr).'


    cdf2.attrs['Author'] = 'Ashley Greeley'
    cdf2.attrs['CreateDate'] = datetime.datetime.now()
    cdf2.close()


if __name__ == '__main__':
    sc_id = 'A'
    day = '2019-01-26'
    #path = load_state(sc_id, day)
