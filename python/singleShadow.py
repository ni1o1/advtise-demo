# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 09:47:02 2022

@author: acer
"""
import pandas as pd
import geopandas as gpd
import os
from suncalc import get_position
from shapely.geometry import Polygon,LineString,Point
import datetime
import math
import numpy as np
#读取shp格式的文件并保存
    



#计算空间直线与平面的交点
def calLinePlaneIntersect(line,plane):
    p1 = line[0],p2 = line[1]
    p1D = plane[0] * p1[0] + plane[1] * p1[1] + plane[2] * p1[2] + plane[3]
    p1D2 = plane[0] * (p2[0] - p1[0]) + plane[1] * (p2[1] - p1[1]) + plane[2] * (p2[2] - p1[2])
    m = abs(p1D / p1D2)
    
    #p2D = (plane[0] * p2[0] + plane[1] * p2[1] + plane[2] * p2[2] + plane[3])
    x = p1[0] - m * (p2[0] - p1[0]);
    y = p1[1] - m * (p2[1] - p1[1])
    z = p1[2] - m * (p2[2] - p1[2])
        
    #if (isNaN(z)):
            #return null
    
    return [x, y, z]

def radianToAngle(radian):
    radian = radian*180/math.pi
    #if radian<0:
     #   radian += 360
    return radian 

def angleToRadian(angle):
    angle = angle/180*math.pi
    #if radian<0:
     #   radian += 360
    return angle 

def lineCrossMultiply(l1, l2, option = "lineSegment") :
    if (option == "lineSegment") :
        #根据两直线计算：前-后
        vectorL1 = [l1[0][0] - l1[1][0], l1[0][1] - l1[1][1], l1[0][2] - l1[1][2]]
        vectorL2 = [l2[0][0] - l2[1][0], l2[0][1] - l2[1][1], l2[0][2] - l2[1][2]]
               
        #计算叉乘
        A = vectorL1[1] * vectorL2[2] - vectorL2[1] * vectorL1[2]
        B = vectorL2[0] * vectorL1[2] - vectorL1[0] * vectorL2[2]
        C = vectorL1[0] * vectorL2[1] - vectorL2[0] * vectorL1[1]

    elif (option == "vector") :
        A = l1[1] * l2[2] - l2[1] * l1[2]
        B = l2[0] * l1[2] - l1[0] * l2[2]
        C = l1[0] * l2[1] - l2[0] * l1[1]
            
    return [A, B, C, A + B + C]

def calPlaneByTwoVectors(n1,n2,p):
    n = lineCrossMultiply(n1, n2, "vector");
    n.splice(3, 1)
    D = -n[0] * p[0] - n[1] * p[1] - n[2] * p[2]
    n.push(D)
    return n


        

#计算每个面的阴影
#shape的格式：两个平面点组成的列表
def calSunShadow(shape,shapeHeight,sunPosition):
    azimuth = (sunPosition['azimuth'])
    altitude = (sunPosition['altitude'])
    symbol = [1,-1]#经度，纬度的符号
   
    if azimuth < 0:
        azimuth += math.pi
        symbol[0] *=-1
    
    distance = shapeHeight/math.tan(altitude)
    
    #计算投影位置偏移
    lonDistance = symbol[0]*distance*math.sin(azimuth)
    latDistance = symbol[1]*distance*math.cos(azimuth)
    
    shadowShape = []
    for i in range(0,2):
        vertex = shape[i]
        #vertex.append(0)
        shadowShape.append(vertex)
        
    for i in range(2,4):    #计算建筑物的顶部点投影位置
        vertex = shape[3-i]
        shadowVertexLon = vertex[0] + lonDistance#经度
        shadowVertexLat = vertex[1] + latDistance#纬度
        shadowShape.append([shadowVertexLon,shadowVertexLat])
    vertex = shadowShape[0]
    shadowShape.append(vertex)
    
    return shadowShape


#多维数据类型：numpy
#输入的shape是一个矩阵（n*2*2) n个建筑物面，每个建筑有2个点，每个点有三个维度
#shapeHeight(n) 每一栋建筑的高度都是一样的
def calSunShadow1(shape,shapeHeight,sunPosition):
    azimuth = (sunPosition['azimuth'])
    altitude = (sunPosition['altitude'])
    symbol = [1,-1]#经度，纬度的符号
   
    if azimuth < 0:
        azimuth += math.pi
        symbol[0] *=-1
        
    distance = shapeHeight/math.tan(altitude)
    
    #计算投影位置偏移
    lonDistance = symbol[0]*distance*math.sin(azimuth)#n个偏移量[n]
    latDistance = symbol[1]*distance*math.cos(azimuth)
    
    n = np.shape(shape)[0]
    shadowShape = np.zeros(n,5,2)#n个建筑物面，每个面都有5个点，每个点都有个维数
   
    shadowShape[:,0:1,:] += shape  #前两个点不变
    shadowShape[:,2:3,0] += shape+ lonDistance 
    shadowShape[:,2:3,1] += shape+ latDistance 
    temp = shadowShape[:,3,:]
    shadowShape[:,3,:] = shadowShape[:,2,:]
    shadowShape[:,2,:] = temp
    
    shadowShape[:,4,:] = shadowShape[:,0,:]
    
    return shadowShape[:,0:1,:]
        




#计算广告牌每个面的阴影
def calBdShadow(shape,shapeHeight,bdPosition,visualArea):
    
    visualGroundR = visualArea['visualGroundR']
    bdHeight = bdPosition[2]
    groundPlane = [0,0,1,0]
    shadowShape = []
    for i in range(0,2):
        vertex = shape[i]
        #vertex.append(0)
        shadowShape.append(vertex)
        
    for i in range(2,3):    #计算建筑物的顶部点投影位置
        
        #如果太阳高度比较高
        if bdHeight < shapeHeight:
            bdHeight = shapeHeight+0.1
        vertex = shape[i]
        vertex.append(shapeHeight)
        shadowVertex = calLinePlaneIntersect([bdHeight,shape[i]],groundPlane)  #计算投影点
        shadowShape.append(shadowVertex[0:1])
    vertex = shadowShape[0]
    shadowShape.append(vertex)

    return  shadowVertex




#def calLineByAngle(azimuth altitude,point):
    
            
if __name__ == "__main__":       
    #data = gpd.read_file(r"G:\\myself\\g\\data\\suzhou\\苏州(1).shp")
    #print(type(data['geometry']))



    #data.set_crs(epsg=3857)
    #data=data['geometry'].to_crs({'init':'epsg:2381'})
    
    lon1 = 139.698311
    lat1 = 35.533796
    lon2 = 139.698311
    lat2 = 35.533642
    #shape = [[ 139.698311,35.533796],[139.698311,35.533642]]
    shape = [[1681062.796660666, 3529393.4372565956],[1681068.8847339682, 3529385.79911048]]
    height = 15
    #给定时间
    #date = datetime.datetime(2022,3,31,0,00,00)
    date = pd.to_datetime('2022-3-3 00:00:00')
    lon = 139.698802
    lat = 35.533816
    sunPosition = get_position(date,lon,lat)#azimuth altitude
    #sunPosition1 = get_position(datetime.datetime(2022,3,31,18,00,00),lon,lat)#azimuth altitude
    #奇怪的现象，11点就过了最南边了
    
    a = calSunShadow(shape,height,sunPosition)
   # b = calSunShadow1(shape,height,sunPosition)
    
    #x = 1
    #y = math.tan(azimuth)*1
    #r = math.sqrt(x*x+y*y)
    #z = math.tan(altitude)*abs(r)
    os.system("pause")







