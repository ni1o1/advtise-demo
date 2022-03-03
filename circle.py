# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 09:20:53 2022

@author: acer
"""

import shapefile
import matplotlib.pyplot as plt
import math
import os
import numpy as np
from osgeo import ogr
from osgeo import osr

'''
0、读取shp数据
1、给定两个点（或者中心点）及其分辨率
2、求该中心点的可视范围及其与地面形成的圆形（shp）
3、将其保存为一个含有地理坐标的shp
4、求建筑物在地面上的投影，赋予不同颜色的值
5、可视区域减去投影区域的叠加
'''
#def get_center_point():
    
def read_shp(sfname):
    '''用shapefile读取shp
    sf = shapefile.Reader(sfname)
    sf_type = sf.shapeType
    sf_bbox = sf.bbox
    sf_numRecords = sf.numRecords
    
    sf_border = sf.shapes()    #提取边界信息
       
    '''
    shp_test = ogr.Open(sfname)
    ''' 也可用下面代码替换，利用driver读取
    driver = ogr.GetDriverByName('ESRI Shapefile')
    sf = driver.open(sfname,0)  #0是只读，1是可写
    if sf is None:
        print('Could not open')
        os.exit(1)
    print(sf)
    '''
    layer = shp_test.GetLayer()
    spatialRef = layer.GetSpatialRef()
    #print(spatialRef)
    #print(dir(layer) ) #获取Layr可以使用的方法
    n = layer.GetFeatureCount()  #获取layer中的元素个数
    print('feature count:',n)
    
    #提取包围图形
    extent = layer.GetExtent()
    feature = layer.GetNextFeature()
    while feature:
        feature = layer.GetNextFeature()
        try:
            area = feature.GetField('shape_area')
            print(area)
        except:
            #print('Done!')
            pass
    layer.ResetReading() #复位
    
    #提取集合要素信息
    feature = layer.GetFeature(1)
    geom = feature.GetGeometryRef()
    dir(geom)
    #关闭文件
    feature.Destroy()
    shp_test.Destroy()
    return spatialRef
        
    #ogr.approximateArcAngles
    #border_points = border.points   #
    
def write_shape(visualArea,sfname,shape,spatialRef):
    '''
    Parameters
    ----------
    visualArea : TYPE
        DESCRIPTION.
    sfname : string
        DESCRIPTION.
    shape : string
        The shape of shp which will be created.

    Returns
    -------
    None.
    '''

    if shape == 'point': #如果需要创建得形状为点
        brandCenterPoint=visualArea['brandCenterPoint']
        height = visualArea['visualHeight']
        
        driver = ogr.GetDriverByName('ESRI Shapefile')
        if os.access(sfname,os.F_OK):   #如果有，删掉
            driver.DeleteDataSource(sfname)
        ds = driver.CreateDataSource(sfname)
        
        #定义空间参考
        #spatialref = osr.SpatialReference()
        #spatialref.ImportFromEPSG(4326)
        
        #创建层
        layer = ds.CreateLayer('point',srs = spatialRef,geom_type = ogr.wkbPoint)  #wkb:文本标记语言二进制形式
        wkt = 'POINT('+ str(brandCenterPoint[0])+' ' +str(brandCenterPoint[1])+')'
    
        #创建属性
        fieldID = ogr.FieldDefn("Height",ogr.OFTString)  #生成一个字符串类型的字符
        fieldID.SetWidth(20)#设置属性显示宽度
        layer.CreateField(fieldID,1)  #
        
        #创建几何图形要素
        feature = ogr.Feature(layer.GetLayerDefn())  #创建要素
        feature.SetField('Height',height)  #设置字段值
        point = ogr.CreateGeometryFromWkt(wkt)  #创建几何要素
        feature.SetGeometry(point)  #使用点
        layer.CreateFeature(feature)   #添加要素到层

        ds.Destroy()   #销毁文件
        
        
        
    if shape == "circle":
        
        #数据
        visualGroundR = visualArea['visualGroundR']
        visualCenter = visualArea['visualCenter']
        visualHeight = visualArea['visualHeight']
        
        driver = ogr.GetDriverByName('ESRI Shapefile')
        if os.access(sfname,os.F_OK):   #如果有，删掉
            driver.DeleteDataSource(sfname)
        ds = driver.CreateDataSource(sfname)
        
        #先创建一个ring    
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring = ogr.ApproximateArcAngles(visualCenter[0],visualCenter[1], 0, 
                                        visualGroundR,visualGroundR, 0.0, 0.0, 360.0, 1.0 )
        ring.CloseRings()
        
        layer = ds.CreateLayer('Polygon',srs = spatialRef,geom_type = ogr.wkbPolygon)
        wkt = 'POLYGON('+ str(ring)[10:]+')'
            
        #创建属性
        fieldID = ogr.FieldDefn("Height",ogr.OFTString)  #生成一个字符串类型的字符
        fieldID.SetWidth(20)#设置属性显示宽度
        layer.CreateField(fieldID,1)  #
        
        #创建几何图形要素
        feature = ogr.Feature(layer.GetLayerDefn())  #创建要素
        feature.SetField('Height',0)  #设置字段值
        point = ogr.CreateGeometryFromWkt(wkt)#生成实体点
        feature.SetGeometry(point)  #使用点
        layer.CreateFeature(feature)   #添加要素到层
        
        ds.Destroy()   #销毁文件

# draw_area(visualArea):
    
def angle2radian(angle):
    return angle*math.pi/180

def radian2angle(radian):
    return radian*180/math.pi    

def cal_visual_area(brandCenterPoint,xResolution = 0.01,isAngle = True,
                eyeResolution = 3,positionAngle=0.0):
    '''

    Parameters
    ----------
    brandCenterPoint : point
        DESCRIPTION.
    xResolution : float, optional
        DESCRIPTION. The default is 0.01.
    isAngle : bool, optional
        angle or radian. The default is True.
    eyeResolution : float, optional
        resolution of eyes("). The default is 3.
    positionAngle : float, optional
        Angle with the real coordinate system(0-360). The default is 0.0.

    Returns
    -------
    visualArea : TYPE
        DESCRIPTION.

    '''
    x,y,z = brandCenterPoint
    if isAngle == True:
        eyeResolution = (eyeResolution/60)/60
        eyeResolution = (eyeResolution*math.pi)/180    #人眼分辨率，弧度
    D = xResolution/eyeResolution
    #半径
    visualR = D/2
    visualGroundR = np.sqrt((D**2)/4 -(z**2))   #地面上的可视化半径
    
    #根据坐标轴角度及半径计算实际位置
    positionAngle = angle2radian(positionAngle) #三角函数计算为弧度制
    
    if positionAngle < 270 and positionAngle >= 90:
        visualCenter = (x-np.cos(positionAngle)*visualGroundR,
                        y-np.sin(positionAngle)*visualGroundR)  #中心点
    elif positionAngle <90 or positionAngle>=270:
       visualCenter = (x+np.cos(positionAngle)*visualGroundR,
                       y+np.sin(positionAngle)*visualGroundR)  #中心点
    visualHeight = z
    
    visualArea = {'brandCenterPoint':brandCenterPoint,
                  'visualR':visualR,
                  'visualGroundR':visualGroundR,
                  'visualCenter':visualCenter,
                  'visualHeight':visualHeight}
    print(visualArea)
    return visualArea
    

if __name__ == "__main__":
    os.chdir(r'G:\myself\g\bd')
    openSfname ='bd.shp'
    
    #读取房屋shp
    inSpayialRef = read_shp(openSfname)
    
    writePointSfname = '18.shp'
    writeCircleSfname = '20.shp'
    brandCenterPoint =[100,39,10] 
    #计算可视区域
    visualArea = cal_visual_area(brandCenterPoint,xResolution = 0.01,isAngle = True,
                                 eyeResolution = 3,positionAngle=30.0)
    write_shape(visualArea,writePointSfname,'point',inSpayialRef)
    write_shape(visualArea,writeCircleSfname,'circle',inSpayialRef)
    
    
