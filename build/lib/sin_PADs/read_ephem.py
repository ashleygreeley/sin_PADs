import h5py
import numpy as np
import datetime
import glob
from os.path import exists
import subprocess
import os
import time

import sin_PADs
from sin_PADs import config

def readephemdata(sc_id,date,mo,day):

        orbleg=[]
        orblegdt=[]
        orbname=[]


        ephemfileglob=glob.glob(sin_PADs.project_data_dir+'/RBSP/ephemeris/'+sc_id+'/rbsp'+sc_id.lower()+'_def_MagEphem_TS04D_'+str(date)+str(mo)+str(day)+'_v*.h5')


        if len(ephemfileglob)==0:
            os.chdir(sin_PADs.project_data_dir+'/RBSP/ephemeris/'+sc_id)
            subprocess.call(['wget','-r','-np','-nd','-nc','-A','rbsp'+sc_id.lower()+'_def_MagEphem_TS04D_'+str(date)+str(mo)+str(day)+'_v*.h5','https://rbsp-ect.newmexicoconsortium.org/data_pub/rbsp'+str(sc_id)+'/MagEphem/definitive/'+str(date)+'/'])
        time.sleep(10)
        ephemfileglob=glob.glob(sin_PADs.project_data_dir+'/RBSP/ephemeris/'+sc_id+'/rbsp'+sc_id.lower()+'_def_MagEphem_TS04D_'+str(date)+str(mo)+str(day)+'_v*.h5')
        while len(ephemfileglob)==0:
            time.sleep(10)
            ephemfileglob=glob.glob(sin_PADs.project_data_dir+'/RBSP/ephemeris/'+sc_id+'/rbsp'+sc_id.lower()+'_def_MagEphem_TS04D_'+str(date)+str(mo)+str(day)+'_v*.h5')
		
        ephemfile=h5py.File(ephemfileglob[0],'r')
        ephemfile=ephemfile
        Aptimes=ephemfile['ApogeeTimes']
        Pertimes=ephemfile['PerigeeTimes']

        for il in range(0,len(Aptimes)):
            Apstr=Aptimes[il].decode('UTF-8')
            Apdatetime=datetime.datetime.strptime(Apstr,'%Y-%m-%dT%H:%M:%S.%fZ')
            Apdectime=float(Apdatetime.strftime("%H"))/24.+float(Apdatetime.strftime("%M"))/1400.+float(Apdatetime.strftime("%S"))/86400.
            orbleg.append(Apdatetime.timetuple().tm_yday+Apdectime)
            orblegdt.append(Apdatetime)
            orbname.append(0)
            if ((il == 0 )):
                Apfirst=Apdatetime.timetuple().tm_yday+Apdectime

        for il in range(0,len(Pertimes)):
            Perstr=Pertimes[il].decode('UTF-8')
            Perdatetime=datetime.datetime.strptime(Perstr,'%Y-%m-%dT%H:%M:%S.%fZ')
            Perdectime=float(Perdatetime.strftime("%H"))/24.+float(Perdatetime.strftime("%M"))/1400.+float(Perdatetime.strftime("%S"))/86400.
            orblegdt.append(Perdatetime)
            orbleg.append(Perdatetime.timetuple().tm_yday+Perdectime)
            orbname.append(1)
            if ((il == 0 )):
                Perfirst=Perdatetime.timetuple().tm_yday+Perdectime
        zipped_lists = zip(orbleg, orbname)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        orbleg, orbname = [ list(tuple) for tuple in  tuples]
        orblegdt.sort()
        if Apfirst < Perfirst:
            whichorb='Apogee'
            firsttime=Apfirst
            ibob=np.append([1],orbname)
        if Perfirst < Apfirst:
            whichorb='Perigee'
            firsttime=Perfirst
            ibob=np.append([0],orbname)
        ibob=np.array(ibob,dtype=object)
        ibob[ibob==0]=str('inbound')
        ibob[ibob==1]='outbound'
        return ibob,orbleg
