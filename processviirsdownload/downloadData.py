#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:15:49 2017

@author: mschull
"""

import os
import subprocess
import h5py
import numpy as np
import pandas as pd
import glob
import datetime
import shutil
from osgeo import gdal,osr
import argparse
import urllib2, base64
from bs4 import BeautifulSoup
import requests
import urllib
from pydap.client import open_url
from pydap.cas import urs
from joblib import Parallel, delayed
import time as timer
import ephem
import sqlite3
import imaplib
import email
import smtplib


ncdcURL = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_pgbh/'
ncdcfluxURL = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_flxf/'

realtimeURL = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/cfs/prod/'

def folders(base):
    data_path = os.path.abspath(os.path.join(base,os.pardir,'VIIRS_DATA'))
    if not os.path.exists(data_path):
        os.makedirs(data_path) 
    static_path = os.path.join(base,"STATIC")
    if not os.path.exists(static_path):
        os.makedirs(static_path) 
    CFSR_path = os.path.join(static_path,"CFSR")   
    if not os.path.exists(CFSR_path):
        os.makedirs(CFSR_path)
    out = {'data_path':data_path,'static_path':static_path,'CFSR_path':CFSR_path}
    return out

base = os.getcwd()
Folders = folders(base)
data_path = Folders['data_path']
static_path = Folders['static_path']
CFSR_path = Folders['CFSR_path']

def warp(args):
    """with a def you can easily change your subprocess call"""
    # command construction with binary and options
    options = ['gdalwarp']
    options.extend(args)
    # call gdalwarp 
    subprocess.check_call(options)

def writeArray2Tiff(data,res,UL,inProjection,outfile,outFormat):

    xres = res[0]
    yres = res[1]

    ysize = data.shape[0]
    xsize = data.shape[1]

    ulx = UL[0] #- (xres / 2.)
    uly = UL[1]# - (yres / 2.)
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outfile, xsize, ysize, 1, outFormat)
    #ds = driver.Create(outfile, xsize, ysize, 1, gdal.GDT_Int16)
    
    srs = osr.SpatialReference()
    
    if isinstance(inProjection, basestring):        
        srs.ImportFromProj4(inProjection)
    else:
        srs.ImportFromEPSG(inProjection)
        
    ds.SetProjection(srs.ExportToWkt())
    
    gt = [ulx, xres, 0, uly, 0, -yres ]
    ds.SetGeoTransform(gt)
    
    ds.GetRasterBand(1).WriteArray(data)
    #ds = None
    ds.FlushCache() 

def is_odd(num):
   return num % 2 != 0

def getGrabTime(time):    
    return int(((time/600)+1)*600)


def getGrabTimeInv(grab_time,doy):
    if is_odd(grab_time):
        hr = grab_time-3
        forecastHR = 3
    else:
        hr = grab_time
        forecastHR = 0
    if hr == 24:
        hr = 0
        doy+=1
    return hr, forecastHR,doy 

def listFD(url, ext=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

def listOrderDir(url, orderID=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [node.get('href') for node in soup.find_all('a') if node.get('href').startswith(orderID)]

def tile2latlon(tile):
    row = tile/24
    col = tile-(row*24)
    # find lower left corner
    lat= (75.-row*15.)-15.
    lon=(col*15.-180.)-15. 
    return [lat,lon]
       

def get_VIIRS_bounds(fn):
    f = h5py.File(fn, 'r')
    east = []
    west = []
    north = []
    south = []
    for i in range(4):
        east = np.append(east,f['Data_Products']['VIIRS-I5-SDR']['VIIRS-I5-SDR_Gran_%d' % i ].attrs['East_Bounding_Coordinate'][0][0])
        west = np.append(west,f['Data_Products']['VIIRS-I5-SDR']['VIIRS-I5-SDR_Gran_%d' % i ].attrs['West_Bounding_Coordinate'][0][0])
        north = np.append(north,f['Data_Products']['VIIRS-I5-SDR']['VIIRS-I5-SDR_Gran_%d' % i ].attrs['North_Bounding_Coordinate'][0][0])
        south = np.append(south,f['Data_Products']['VIIRS-I5-SDR']['VIIRS-I5-SDR_Gran_%d' % i ].attrs['South_Bounding_Coordinate'][0][0])
    east = east.min()
    west = west.max()
    north = north.max()
    south = south.min()
    date = f['Data_Products']['VIIRS-I5-SDR']['VIIRS-I5-SDR_Gran_0'].attrs['Beginning_Date'][0][0]
    year = int(date[:4])
    month = int(date[4:6])
    day = int(date[6:])
    dd = datetime.datetime(year,month,day)
    doy = ((dd-datetime.datetime(year,1,1)).days)+1
    N_Day_Night_Flag = f['Data_Products']['VIIRS-I5-SDR']['VIIRS-I5-SDR_Gran_0'].attrs['N_Day_Night_Flag'][0][0]
    bounds1 = {'filename':[fn],'N_Day_Night_Flag':[N_Day_Night_Flag]}
    bounds2 = {'east':[east],'west':[west],'north':[north],'south':[south]}
    bounds3 = {'doy':[doy],'year':[year]}

    df1 = pd.DataFrame.from_dict(bounds1) 
    df2 = pd.DataFrame.from_dict(bounds2) 
    df3 = pd.DataFrame.from_dict(bounds3)
    dictDF = pd.concat([df1,df2,df3],axis=1,copy=False)
    return dictDF

def downloadSubscriptionSDR(inurl=None): 
#    if year==None:
#        dd = datetime.date.today()+datetime.timedelta(days=-1)
#        year = dd.year
#    if doy==None:
#        dd = datetime.date.today()+datetime.timedelta(days=-1)
#        month = dd.month
#        day = dd.day
#    else:
#        dd=datetime.datetime(year,1,1)+datetime.timedelta(days=doy)
#        month = dd.month
#        day = dd.day
#    filePath = os.path.join(data_path,"%d" % year,"%02d" % month)
    #download I5 data
    ext = 'h5'
    if inurl==None: # Use subscription
        url = 'https://download.class.ngdc.noaa.gov/download/sub/hain/85113/'
        years = []
        months = []
        days = []
        for fn in listFD(url, ext):
            fileName = str(fn.split('/')[-1])  
#            if (fileName.split("_")[2]=='d%d%02d%02d' % (year,month,day)):
            year = int(fileName.split("_")[-7][1:5])
            month = int(fileName.split("_")[-7][5:7])
            day = int(fileName.split("_")[-7][7:9])
            years.append(year)
            months.append(month)
            days.append(day)
            filePath = os.path.join(data_path,"%d" % year,"%02d" % month)
            if not os.path.exists(filePath):
                os.makedirs(filePath)
        
            outName=os.path.join(filePath,fileName)
            
            if not os.path.isfile(outName):
                print "downloading:  %s" % fileName
                urllib.urlretrieve(url+fileName, outName)

    # download cloud data
    if inurl==None:
        url = 'https://download.class.ngdc.noaa.gov/download/sub/hain/85123/'
        years = []
        months = []
        days = []
        for fn in listFD(url, ext):
            fileName = str(fn.split('/')[-1])  
#            if (fileName.split("_")[2]=='d%d%02d%02d' % (year,month,day)):
            year = int(fileName.split("_")[-7][1:5])
            month = int(fileName.split("_")[-7][5:7])
            day = int(fileName.split("_")[-7][7:9])
            years.append(year)
            months.append(month)
            days.append(day)
            filePath = os.path.join(data_path,"%d" % year,"%02d" % month)
            if not os.path.exists(filePath):
                os.makedirs(filePath)
        
            outName=os.path.join(filePath,fileName)
            
            if not os.path.isfile(outName):
                print "downloading:  %s" % fileName
                urllib.urlretrieve(url+fileName, outName)
    else:
        years = []
        months = []
        days = [] 
        for fn in listFD(inurl, ext):
            fileName = str(fn.split('/')[-1])  
#            if (fileName.split("_")[2]=='d%d%02d%02d' % (year,month,day)):
            year = int(fileName.split("_")[-7][1:5])
            month = int(fileName.split("_")[-7][5:7])
            day = int(fileName.split("_")[-7][7:9])
            years.append(year)
            months.append(month)
            days.append(day)
            filePath = os.path.join(data_path,"%d" % year,"%02d" % month)
            if not os.path.exists(filePath):
                os.makedirs(filePath)
        
            outName=os.path.join(filePath,fileName)
            
            if not os.path.isfile(outName):
                print inurl+fileName
                print "downloading:  %s" % fileName
                urllib.urlretrieve(inurl+fileName, outName)
                
    date_df = pd.DataFrame.from_dict({'years': years,'months': months,'days': days})
    
    return date_df                
                
                
def getInsolation(earthLoginUser,earthLoginPass,tile,year=None,doy=None):
    if year==None:
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        year = dd.year
    if doy==None:
        doy = (datetime.date.today()-datetime.date(year,1,1)).days
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        month = dd.month
        day = dd.day
    llLat,llLon = tile2latlon(tile)
    ulx = llLon
    uly = llLat+15.
    lrx = llLon+15.
    lry = llLat
    MERRA2_ulLat = 90.0
    MERRA2_ulLon = -180.0
    MERRA2LatRes = 0.5
    MERRA2LonRes = 0.625
    nrow = 3750
    ncol = 3750
    date = '%d%03d' % (year,doy)
    inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    #========get MERRA2 Insolation data at overpass time====================== 
    inRes = [MERRA2LonRes,MERRA2LatRes]
    inUL = [MERRA2_ulLon,MERRA2_ulLat]

    if year <1992:
        fileType = 100
    elif year >1991 and year < 2001:
        fileType=200
    elif year > 2000 and year<2011:
        fileType = 300
    else:
        fileType = 400

    opendap_url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/'
    product = 'M2T1NXRAD.5.12.4'
    
    filename = 'MERRA2_%d.tavg1_2d_rad_Nx.%04d%02d%02d.nc4' % (fileType,year,month,day)
    fullUrl =os.path.join(opendap_url,product,'%04d'% year,'%02d'% month,filename)
    session = urs.setup_session(username = earthLoginUser, 
                password = earthLoginPass,
                check_url=fullUrl)
    d = open_url(fullUrl,session=session)
    Insol = d.SWGDNCLR
        
    #====get daily insolation=========================================
    outFN = os.path.join(static_path,'INSOL24','T%03d' % tile, 'RS24_%s_T%03d.tif' % (date,tile))
    if not os.path.exists(outFN):
        dataset2 = np.flipud(np.sum(np.squeeze(Insol[:,:,:]),axis=0))
        outfile = os.path.join(os.getcwd(),'insol24')
        outFormat = gdal.GDT_Float32
        writeArray2Tiff(dataset2,inRes,inUL,inProj4,outfile,outFormat)
        optionList = ['-overwrite', '-s_srs', '%s' % inProj4,'-t_srs','%s' % inProj4,\
        '-te', '%f' % ulx, '%f' % lry,'%f' % lrx,'%f' % uly,'-r', 'bilinear',\
        '-ts', '%f' % nrow, '%f' % ncol,'-multi','-of','GTiff','%s' % outfile, '%s' % outFN]
        warp(optionList)

def downloadCFSRpython(hr1file,year=None,doy=None):  
    if year==None:
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        year = dd.year
        
    if doy==None:
        doy = (datetime.date.today()-datetime.date(year,1,1)).days
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        month = dd.month
        day = dd.day
        url = realtimeURL+'cdas.%d%02d%02d/' % (year,month,day)
    else:
        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(doy-1)
        if (datetime.date.today()-datetime.date(year,dd.month,dd.day)).days > 7:
            url = os.path.join(ncdcfluxURL,"%s" % year,"%d%02d" % (year,dd.month),
                                "%d%02d%02d" % (year,dd.month,dd.day))
        else:
            url = realtimeURL+'cdas.%d%02d%02d/' % (year,dd.month,dd.day)
    dstpath =  os.path.join(CFSR_path,"%d" % year,"%03d" % doy)
    if not os.path.exists(dstpath):
        os.makedirs(dstpath)  
#    hrs = (HRs/6)*HRs[0]
#    forecastHRs = HRs-hrs
#    for i in range(len(HRs)):
#
##    for i in range(0,8):
##        hrs = [0,0,6,6,12,12,18,18]
#        hr = hrs[i]
##        forcastHRs = [0,3,0,3,0,3,0,3]
#        forecastHR = forecastHRs[i]
#        hr1file = 'cdas1.t%02dz.sfluxgrbf%02d.grib2' % (hr,forecastHR)
    print url
    pydapURL = os.path.join(url,hr1file)
    outFN = os.path.join(dstpath,hr1file)
    if not os.path.exists(outFN):
        print "downloading file...%s" % hr1file
        getHTTPdata(pydapURL,outFN)

def getCFSRInsolation(tile,year=None,doy=None):
     
    if year==None:
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        year = dd.year
        
    if doy==None:
        doy = (datetime.date.today()-datetime.date(year,1,1)).days
        dd = datetime.date.today()+datetime.timedelta(days=-1)
    else:
        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(doy-1)

        
    dstpath =  os.path.join(CFSR_path,"%d" % year,"%03d" % doy)
    if not os.path.exists(dstpath):
        os.makedirs(dstpath)
        
    tile_path = os.path.join(dstpath,"T%03d" % tile)
    if not os.path.exists(tile_path): 
        os.makedirs(tile_path)
    date = '%d%03d' % (year,doy)
    outDIR = os.path.join(static_path,'INSOL24')
    if not os.path.exists(outDIR): 
        os.makedirs(outDIR)
    outFN = os.path.join(outDIR, 'RS24_%s_T%03d.tif' % (date,tile))
    if not os.path.exists(outFN):    # If it exists skip
        llLat,llLon = tile2latlon(tile)
        ulx = llLon
        uly = llLat+15.
        lrx = llLon+15.
        lry = llLat
        nrow = 3750
        ncol = 3750
        inProj4 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    
        print "date:%s" % date
        print("tile:T%03d" % tile)
        print "============================================================"
    
        o = ephem.Observer()
        o.lat, o.long = '%3.2f' % (llLat+7.5), '%3.2f' % (llLon+7.5)
        sun = ephem.Sun()
        #================finding the local noon====================================
        dd = datetime.datetime(dd.year,dd.month,dd.day,0,0)
        sunrise = o.previous_rising(sun, start=dd)
        noon = o.next_transit(sun, start=sunrise)
        hr = noon.datetime().hour
        #================finding the sunrise and sunset ===========================
        dd = datetime.datetime(dd.year,dd.month,dd.day,hr,0)
        sunrise = o.previous_rising(sun, start=dd)
        t_rise = sunrise.datetime().hour
        sunset = o.next_setting(sun, start=dd)
        t_end = sunset.datetime().hour
        doy_end = (sunset.datetime()-datetime.datetime(year,1,1)).days+1
        
        grab_time = getGrabTime((t_rise-1)*100)
        firstHR = getGrabTimeInv(grab_time/100,doy)
        grab_time = getGrabTime((t_end+1)*100)
        lastHR = getGrabTimeInv(grab_time/100,doy_end)
        doy_end = lastHR[2]
        print doy, doy_end
        print firstHR, lastHR
        
        if doy_end>doy:
            HRs=[]
            doys=[]
            d1HRs = np.array(range(firstHR[0],24,3))
            doy1 = np.tile(doy,len(d1HRs))
            d2HRs = np.array(range(0,lastHR[0],3))
            if len(d2HRs)==0:
                d2HRs=np.array([0])
            doy2 = np.tile(doy_end,len(d2HRs))
            HRs=np.append(d1HRs,d2HRs)
            doys = np.append(doy1,doy2)
        else:
            HRs = np.array(range(firstHR[0],lastHR[0]+1,3))
            doys = np.tile(doy,len(HRs))
        hrs = (HRs/6)*6
        forecastHRs = HRs-hrs
        outData = []
        print HRs
        for i in range(len(HRs)):
    #        hrs = [0,0,6,6,12,12,18,18]
            hr = hrs[i]
            doy = doys[i]
    #        forcastHRs = [0,3,0,3,0,3,0,3]
            tile_path = os.path.join(dstpath,"T%03d" % tile)
            if not os.path.exists(tile_path): 
                os.makedirs(tile_path)
            dstpath =  os.path.join(CFSR_path,"%d" % year,"%03d" % doy)
            if not os.path.exists(dstpath):
                os.makedirs(dstpath)
            forecastHR = forecastHRs[i]
            hr1file = 'cdas1.t%02dz.sfluxgrbf%02d.grib2' % (hr,forecastHR)
            downloadCFSRpython(hr1file,year,doy)
    #        print url
    #        pydapURL = os.path.join(url,hr1file)
            gribFN = os.path.join(dstpath,hr1file)
    
            cfsr_out = os.path.join(tile_path, "CFSR_INSOL_%d%03d_%02d00_00%d.tif" % (year,doy,hr,forecastHR))
            cfsr_outvrt = os.path.join(tile_path,"CFSR_INSOL_%d%03d_%02d00_00%d.vrt" % (year,doy,hr,forecastHR))
    #        if not os.path.exists(outFN):
    #            print "downloading file...%s" % hr1file
    #            getHTTPdata(pydapURL,outFN)
                
                #------download file                
        
    #        print "processing file...%s" % hr1file   
            dataset = gdal.Open(gribFN, gdal.GA_ReadOnly)
            ##=====find correct dataset from grib===============
            for bandNum  in range(90,100): # use this range because file is towrds the end
                band = dataset.GetRasterBand(bandNum)
                if (band.GetMetadata_List()[0]=='GRIB_COMMENT=Downward Short-Wave Rad. Flux [W/(m^2)]'):
                    subprocess.check_output(["gdal_translate", "-of", "vrt", "-b",
                                             "%d" % bandNum,"-r", "bilinear",
                                             "%s" % gribFN, "%s" % cfsr_outvrt])
                    subprocess.check_output(["gdalwarp", "-overwrite", "-t_srs", 
                                             '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
                                             "%s" % cfsr_outvrt, "%s" % cfsr_out,"-wo",
                                             "SOURCE_EXTRA=1000","--config",
                                             "CENTER_LONG", "0","-te", "%d"  % ulx,
                                             "%d" % lry,"%d" % lrx, "%d" % uly,"-tr", 
                                             "0.004", "0.004"])
            
    #            gdalwarp -t_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' test.vrt test3.tif -wo SOURCE_EXTRA=1000 --config CENTER_LONG 0 -te -15.0 30.0 0.0 45.0 -tr 0.004 0.004
            g = gdal.Open(cfsr_out, gdal.GA_ReadOnly)
            data = g.ReadAsArray()
            outData.append(np.reshape(data,[nrow*ncol]))
           
        aa = np.array(outData)
        rs24 = np.sum(aa,axis=0)
        rs24 = np.reshape(rs24,[nrow, ncol])
        outFormat = gdal.GDT_Float32
        inUL = [ulx,uly]
        inRes = [0.004,0.004]
        writeArray2Tiff(rs24,inRes,inUL,inProj4,outFN,outFormat)
        print "finished INSOL!"
    
def write_gen_sfc_prof(path,date,hr):
    fn = os.path.join(path,'gen_sfc_prof_fields.gs')
    file = open(fn, "w")
    file.write("'open current.ctl'\n")
    file.write("\n")
    file.write("'set lon -180 180'\n")
    file.write("'set lat -89.875 89.875' \n")
    file.write("\n")
    file.write("'set gxout fwrite'\n")
    file.write("'set fwrite sfc_temp_%s_%02d00.dat'\n" % (date,hr))
    file.write("'d re(smth9(tmpsig995),0.25,0.25)'\n")
    file.write("'disable fwrite'\n")
    file.write("\n")
    file.write("'set gxout fwrite'\n")
    file.write("'set fwrite sfc_pres_%s_%02d00.dat'\n" % (date,hr))
    file.write("'d re(smth9(pressfc),0.25,0.25)'\n")
    file.write("'disable fwrite'\n")
    file.write("\n")
    file.write("'set gxout fwrite'\n")
    file.write("'set fwrite sfc_spfh_%s_%02d00.dat'\n" % (date,hr))
    file.write("'d re(smth9(spfhhy1),0.25,0.25)'\n")
    file.write("'disable fwrite'\n")
    file.write("\n")
    file.write("'set gxout fwrite'\n")
    file.write("'set fwrite sfc_lwdn_%s_%02d00.dat'\n" % (date,hr))
    file.write("'d re(smth9(dlwrfsfc),0.25,0.25)'\n")
    file.write("'disable fwrite'\n")
    file.write("\n")
    file.write("'set gxout fwrite'\n")
    file.write("'set fwrite temp_profile_%s_%02d00.dat'\n" % (date,hr))
    file.write("z=1\n")
    file.write("while (z <= 21)\n")
    file.write(" 'set z 'z\n")
    file.write(" 'd re(smth9(tmpprs),0.25,0.25)'\n")
    file.write(" z=z+1\n")
    file.write("endwhile\n")
    file.write("'disable fwrite'\n")
    file.write("\n")
    file.write("'set gxout fwrite'\n")
    file.write("'set fwrite spfh_profile_%s_%02d00.dat'\n" % (date,hr))
    file.write("z=1\n")
    file.write("while (z <= 21)\n")
    file.write(" 'set z 'z\n")
    file.write(" 'd re(smth9(spfhprs),0.25,0.25)'\n")
    file.write(" z=z+1\n")
    file.write("endwhile\n")
    file.write("'disable fwrite'\n")
    file.close()


def runGrads(fn,gspath):
    g2ctl = 'g2ctl'
    gribmap = 'gribmap'
    opengrads = 'opengrads'
    gen_sfc_prof = os.path.join(gspath,'gen_sfc_prof_fields.gs')
    subprocess.check_output("%s %s > ./current.ctl" % (g2ctl,fn), shell=True)
    subprocess.check_output("%s -i ./current.ctl -0" % gribmap, shell=True)
    out = subprocess.check_output("%s -blxc 'run %s '" % (opengrads,gen_sfc_prof), shell=True)
    
    print out


class earthDataHTTPRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return urllib2.HTTPRedirectHandler.http_error_302(self, req, fp, code, msg, headers)
    

def getHTTPdata(url,outFN,auth=None):
    request = urllib2.Request(url) 
    if not (auth == None):
        username = auth[0]
        password = auth[1]
        base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')
        request.add_header("Authorization", "Basic %s" % base64string) 
    
    cookieprocessor = urllib2.HTTPCookieProcessor()
    opener = urllib2.build_opener(earthDataHTTPRedirectHandler, cookieprocessor)
    urllib2.install_opener(opener) 
    r = opener.open(request)
    result = r.read()
    
    with open(outFN, 'wb') as f:
        f.write(result)
def moveFiles(sourcepath,destinationpath,date,hr):            
    files = glob.glob(os.path.join(sourcepath,"*_%s_%02d00.dat" % (date,hr)))
    for fn in files:
        infn = os.path.join(sourcepath,fn.split(os.sep)[-1])
        outfn = os.path.join(destinationpath,fn.split(os.sep)[-1])
        if fn.split('_')[1]=='profile':
            prof = np.fromfile(infn, dtype=np.float32)
            prof = np.flipud(prof.reshape([21,720,1440]))
            prof.tofile(outfn)
        else:                
            shutil.move(infn, outfn) 

def getCFSRdata(year=None,doy=None):  
    if year==None:
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        year = dd.year
        
    if doy==None:
        doy = (datetime.date.today()-datetime.date(year,1,1)).days
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        month = dd.month
        day = dd.day
        url = realtimeURL+'cdas.%d%02d%02d/' % (year,month,day)
    else:
        dd = datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1)
        url = os.path.join(ncdcURL,"%s" % year,"%d%02d" % (year,dd.month),
                                "%d%02d%02d" % (year,dd.month,dd.day))
    

    dstpath =  os.path.join(CFSR_path,"%d" % year,"%03d" % doy)
    if not os.path.exists(dstpath):
        os.makedirs(dstpath) 
    levs="(100|150|200|250|300|350|400|450|500|550|600|650|700|750|800|850|900|925|950|975|1000) mb"

    
    s1="(HGT):%s" % levs
    s2="(TMP):%s" % levs
    s3="(SPFH):%s" % levs
    s4="DLWRF:surface"
    s5="HGT:surface"
    s6="PRES:surface"
    s7="SPFH:1 hybrid level"
    s8="TMP:0.995 sigma level"
    s9="UGRD:0.995 sigma level"
    s10="VGRD:0.995 sigma level"
    wgrib = "wgrib2"

#    for year in range(iyear,eyear):
#        for doy in range(iday,eday):
    
    date = "%d%03d" %(year,doy)
    print "date:%s" % date
    print "============================================================"
    for i in range(0,8):
        hrs = [0,0,6,6,12,12,18,18]
        hr = hrs[i]
        write_gen_sfc_prof(dstpath,date,hr)
        forcastHRs = [0,3,0,3,0,3,0,3]
        forcastHR = forcastHRs[i]
        hr1file = 'cdas1.t%02dz.pgrbh%02d.grib2' % (hr,forcastHR)
        
        
        
        #------download file                
        pydapURL = os.path.join(url,hr1file)
        outFN = os.path.join(dstpath,hr1file)
        cfsr_out = os.path.join(dstpath,"CFSR_%d%03d_%02d00_00%d.grib2" % (year,doy,hr,forcastHR))
        if not os.path.exists(cfsr_out):
            print "processing file...%s" % hr1file
            print "downloading...%s" % pydapURL
            getHTTPdata(pydapURL,outFN)
            
            #------extract data
    
    
            
            
            subprocess.check_output(["%s" % wgrib, "%s" % outFN, "-match",
                                     "\"%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\"" % (s1,s2,s3,s4,s5,s6,s7,s8,s9,s10),
                                     "-grib", "%s" % cfsr_out])
            os.remove(outFN)
            #-------process using grads------
            runGrads(cfsr_out,dstpath)    
            
    #        moveFiles(os.getcwd(),dstpath,date,forcastHR,hr,"dat")
            moveFiles(os.getcwd(),dstpath,date,hr)
    print "finished processing!"
def createDB(year=None,doy=None):
    
    if year==None:
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        year = dd.year
    if doy==None:
        dd = datetime.date.today()+datetime.timedelta(days=-1)
        month = dd.month
    else:
        dd=datetime.datetime(year,1,1)+datetime.timedelta(days=doy-1)
        month = dd.month

#        
#    df = pd.DataFrame()
##    parDir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
#    dirpath = os.path.join(data_path,"%d" % year, "%02d" % month)
#    database = os.path.join(dirpath,'I5_database.csv')
#    fileList = glob.glob(os.path.join(dirpath, "*SVI05*.h5"))
#    db = pd.read_csv( database )
#    
#    for fn in fileList:
#        filename = fn.split(os.sep)[-1]
#        if (np.sum(db.filename==filename)==0):            
#            try:
#                df1 = get_VIIRS_bounds(os.path.join(dirpath, filename))
#                df = df.append(df1, ignore_index=True)
#            except: 
#              pass
#            
#    for dirpath, dirnames, filenames in os.walk(parDir):
#        try:
#            for filename in [f for f in filenames if ((f.split("_")[0].split(".")[1] == "SVI05") and f.endswith(".h5"))]:
#                try:
#                    df1 = get_VIIRS_bounds(os.path.join(dirpath, filename))
#                    df = df.append(df1, ignore_index=True)
#                except: 
#                  pass
#        except:
#            for filename in [f for f in filenames if (f.startswith("SVI05") and f.endswith(".h5"))]:
#                try:
#                    df1 = get_VIIRS_bounds(os.path.join(dirpath, filename))
#                    df = df.append(df1, ignore_index=True)
#                except: 
#                  pass
    
#    df.to_csv(database, index=False)
    
    #====use sqlite db=========
    I5_db_name = os.path.join(data_path,'viirs_database.db')
    dirpath = os.path.join(data_path,"%d" % year, "%02d" % month)
    fileList = glob.glob(os.path.join(dirpath, "*SVI05*.h5"))
    if not os.path.exists(I5_db_name):
        i5_df = pd.DataFrame()
        conn = sqlite3.connect( I5_db_name )
    else:    
        conn = sqlite3.connect( I5_db_name )
        i5_df = pd.read_sql_query("SELECT * from i5",conn)
    df = pd.DataFrame()
    for fn in fileList:
        filename = fn.split(os.sep)[-1]          
        try:
            df1 = get_VIIRS_bounds(os.path.join(dirpath, filename))
            df = df.append(df1, ignore_index=True)
        except: 
          pass
    i5_df = i5_df.append(df,ignore_index=True)
    i5_df = i5_df.drop_duplicates(subset='east',keep='last')
    i5_df.to_sql("i5", conn, if_exists="replace", index=False)
    
    conn.close()
    print "all done!!"
    
start = timer.time()

def runProcess(tiles,downloadurl=None):
#    if year==None: # if None assume its real-time processing 
#        downloadSubscriptionSDR()
#        getCFSRdata() 
#        for tile in tiles:
#            getCFSRInsolation(tile)
    if downloadurl==None:
        date_df = downloadSubscriptionSDR()
        for i in range(len(date_df)):
            year = date_df['years'][i]
            month = date_df['months'][i]
            day = date_df['days'][i]
            ss = datetime.date(year,month,day)-datetime.date(year,1,1)
            doy = ss.days+1
            getCFSRdata(year,doy)
            for tile in tiles:
                getCFSRInsolation(tile,year,doy)
    else:
        date_df = downloadSubscriptionSDR(downloadurl)
        for i in range(len(date_df)):
            year = date_df['years'][i]
            month = date_df['months'][i]
            day = date_df['days'][i]
            ss = datetime.date(year,month,day)-datetime.date(year,1,1)
            doy = ss.days+1
            getCFSRdata(year,doy)
            for tile in tiles:
                getCFSRInsolation(tile,year,doy)
 
def read_email_from_gmail(emailadd,password):
    try:
        mail = imaplib.IMAP4_SSL("imap.gmail.com")
        mail.login(emailadd,password)
        mail.select('inbox')

#        type, data = mail.search(None, 'ALL')
        type, data = mail.search(None, '(UNSEEN)')
        mail_ids = data[0]
        mail_ids = map(int, mail_ids.split(" "))
        if mail_ids == '':
            classOrderIDs = []
            url = ''
        else:

#        id_list = mail_ids.split()   
#        first_email_id = int(id_list[0])
#        latest_email_id = int(id_list[-1])


            for mail_id in mail_ids:
                typ, data = mail.fetch(mail_id, '(RFC822)' )
                classOrderIDs = []
                url = ''
                for response_part in data:
                    if isinstance(response_part, tuple):
                        msg = email.message_from_string(response_part[1])
                        email_subject = msg['subject']
                        email_from = msg['from']
                        if email_subject.split()[-1] == "Complete":                        
                            classOrderIDs.append(email_subject.split()[-3])
                            print 'From : ' + email_from + '\n'
                            print 'Subject : ' + email_subject + '\n'
                            msg.get_payload()
                            a = msg.get_payload()
                            b = a.splitlines()
                            bb = np.char.strip(b)
                            loc = np.argwhere(np.array(b) == "Alternatively, you can also pick up your data  via")+2
                            url = bb[loc[0][0]]
                        else:
                            print("nothing to see here")

    except Exception, e:
        print str(e)
    return classOrderIDs,url

def sendEmail(orderID):
    fromaddr = "ordersatdata@gmail.com"
    toaddrs = "bucricket@gmail.com"
    #fromaddr = 'user_me@gmail.com'
    #toaddrs  = 'user_you@gmail.com'
    msg = "\r\n".join([
      "From: %s" % fromaddr,
      "To: %s" % toaddrs,
      "Subject: Finished processing order %s" % orderID,
      "",
      ""
      ])
    username = 'ordersatdata@gmail.com'
    password = 'sushmaMITCH12'
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.starttls()
    server.login(username,password)
    server.sendmail(fromaddr, toaddrs, msg)
    server.quit()

def main():
    # Get time and location from user
    parser = argparse.ArgumentParser()
#    parser.add_argument("year", nargs='?', type=int, default=None, help="year of data")
#    parser.add_argument("start_doy", nargs='?',type=int, default=None, help="start day of processing. *Note: leave blank for Real-time")
#    parser.add_argument("end_doy", nargs='?',type=int, default=None, help="end day of processing. *Note: leave blank for Real-time")
#    parser.add_argument('-p','--parentDir', nargs='*',type=int, default=None, help="parent director for large orders from CLASS e-mail. *Note: leave blank for Real-time")
    parser.add_argument('-o','--orderIDs', nargs='*',type=str, default=None, help="list of order IDs from CLASS e-mail. *Note: leave blank for Real-time")
    parser.add_argument('-t','--tiles', nargs='*',type=int, default=None, help='list of tiles')
    parser.add_argument('-c','--cron', nargs='*',type=int, default=0, help='1 for using crontab')
    args = parser.parse_args()
    orderIDs = args.orderIDs
    cron = args.cron
    if cron[0]==1:
        print("crontab")
        orderIDs,url = read_email_from_gmail("ordersatdata@gmail.com","sushmaMITCH12")
#    parentDir = args.parentDir
    tiles = args.tiles
#    if start_doy == None:
#        tiles = [60,61,62,63,64,83,84,85,86,87,88,107,108,109,110,111,112]
#        start = timer.time()
#        runProcess(tiles)
#        createDB()
#        end = timer.time()
#        print("program duration: %f minutes" % ((end - start)/60.))
    if orderIDs ==None: # Use subscription
        if tiles==None:
            tiles = [60,61,62,63,64,83,84,85,86,87,88,107,108,109,110,111,112]
        start = timer.time()
        runProcess(tiles)
        createDB()
        end = timer.time()
        print("program duration: %f minutes" % ((end - start)/60.)) 
    else:
        if tiles==None:
            tiles = [60,61,62,63,64,83,84,85,86,87,88,107,108,109,110,111,112]
        for orderID in orderIDs:
#            url = 'https://download.class.ngdc.noaa.gov/download/%s/' % orderID
#            for order in listOrderDir(url, orderID):
                #download_url = 'https://download.class.ngdc.noaa.gov/download/%s/' % str(order)
            download_url = url+"/"
            if not download_url.split("/")[-2] == '001':
                download_url = download_url+"001/"
            print download_url
            start = timer.time()
            runProcess(tiles,download_url)
            end = timer.time()
            createDB()
            if not orderID == '':
                sendEmail(orderID)
            
            
            print("program duration: %f minutes" % ((end - start)/60.))
#            else:
#                download_url = 'https://download.class.ngdc.noaa.gov/download/%s' % orderID
                
#            if parentDir==None:
##                download_url = 'https://download.class.ncdc.noaa.gov/download/%d/001/' % orderID
##                if not listFD(download_url, 'h5'):
#                download_url = 'https://download.class.ngdc.noaa.gov/download/%d/001' % orderID
#            else:
#                parentDir = args.parentDir[0]
##                download_url = 'https://download.class.ncdc.noaa.gov/download/%d/%d/001/' % (parentDir,orderID)
##                if not listFD(download_url, 'h5'):
#                download_url = 'https://download.class.ngdc.noaa.gov/download/%d/%d/001' % (parentDir,orderID)

#                days = range(start_doy,end_doy)
#                print download_url
#                start = timer.time()
#                for doy in days:
#                    runProcess(tiles,year,doy,download_url)
#                createDB()
#                end = timer.time()
#                print("program duration: %f minutes" % ((end - start)/60.))       
   
main()      


