function angleToRadian(angle) {
    return angle * Math.PI / 180;
};

function radianToAngle(radian) {
    return radian * 180 / Math.PI;
};

//计算可视区域形状
function calObservedArea(observerCenter, observerHeight, observedR) {

    var observedGroundR = Math.sqrt((Math.pow(observedR, 2)) - (Math.pow(observerHeight, 2))); //地面上的可视化半径

    //坐标转换
    //var visualHeight = z;
    observedR /= 1000;
    observedGroundR /= 1000;
    //observerHeight /= 1000;
    observerCenter = turf.toWgs84(observerCenter);

    var observedArea = {
        'observerCenter': observerCenter, //观察者中心位置
        'observerHeight': observerHeight, //观察者高度
        'observedR': observedR, //观察到的半径
        'observedGroundR': observedGroundR, //地面半径
    }
    //console.log("observedArea = ", observedArea);
    return observedArea;
}

//获取地面的可视坐标点
function getCirclePosition(observedArea) {

    //console.log(visualArea);
    var options = { steps: 360, units: 'kilometers' };

    //方法一：生成圆弧线转面
    //var circleLine = turf.lineArc(visualArea.visualCenter, visualArea.visualGroundR, 0, 360, options); //计算圆所有的点
    //console.log("circleLine = ",circleLine);
    //var circlePoly = turf.polygonize(circleLine); //线段转多边形
    //console.log(circleLine,"circlePoly = ",circlePoly);

    //生成圆弧面sector 
    var circlePoly = turf.circle(observedArea.observerCenter, observedArea.observedGroundR, options);
    //console.log("circlePoly = ", circlePoly);

    return circlePoly;
}
//Step1:
//获取可视范围内的建筑
function getVisualBuilding(circlePoly, buildings) {
    //var total = 0;
    var visualBuilding = [];
    var idSelected = [];
    //var exceptBuildingsPoly = circlePoly;

    var start = new Date().getTime();
    //统计有多少个多边形相交
    turf.geomEach(buildings, function (currentGeometry, featureIndex, featureProperties) {
        if (turf.booleanOverlap(circlePoly, currentGeometry) || turf.booleanWithin(currentGeometry, circlePoly)) { //如果包含或者相交

            //total++; //统计个数
            currentGeometry.properties = featureProperties; //特性移植
            //console.log(currentGeometry);

            visualBuilding.push(currentGeometry); //记录可视范围内的建筑
            idSelected.push(currentGeometry.properties.id); //记录可视范围内的建筑编号

        };
    });

    var multiBuildPoly = turf.multiPolygon(visualBuilding); //多个多边形

    var visualBuildings = {
        "multiBuildPoly": multiBuildPoly, //视野内的建筑
        "idSelected": idSelected, //视野内的建筑ID
    };

    var end = new Date().getTime();
    //console.log('cost is', `${end - start}ms`);
    return visualBuildings;
}

//构造函数
buildingsShapes = function (shape, height, angle, shelter, shelterShape) {
    this.shape = shape; //没有特定顺序
    this.height = height;
    this.angle = angle;
    this.shelter = shelter; //遮挡它的面
    this.shelterShape = shelterShape;
}
buildingsShapes.prototype = {
    maxDistance: function (point) {
        //console.log('aka',turf.distance(point, this.shape[0]), turf.distance(point, this.shape[1]));
        var maxDistance = Math.max(turf.distance(point, this.shape[0]), turf.distance(point, this.shape[1]));
        //console.log(maxDistance);
        return maxDistance;
    },
    minDistance: function (point) {
        var minDistance = Math.min(turf.distance(point, this.shape[0]), turf.distance(point, this.shape[1]));
        //console.log(maxDistance);
        return minDistance;
    },
    length: function () {
        return turf.distance(this.shape[0], this.shape[1]);
    },
    isShelter: function (bS, point) { //被这个面遮挡
        //console.log(this.distance(point),bS.distance(point));
        if ((bS.angle[0] < this.angle[1]) && (bS.angle[1] > this.angle[0]) && (this.maxDistance(point) > bS.minDistance(point))) {
            this.shelter.push(bS);
            return true;
        } else { return false; }
    },
}

//Step2:
//获取建筑物的每个面的三维坐标集合
//计算可视范围内的建筑物面
function calVisibleBuilding(visualBuildings, observedArea) { //circlePoly, 
    var buildings = visualBuildings.multiBuildPoly.geometry.coordinates;
    var observerCenter = observedArea.observerCenter;
    var shapesGroup = [];

    //console.log(buildings);
    for (var i = 0; i < buildings.length; i++) {
        var currentBuilding = visualBuildings.multiBuildPoly.geometry.coordinates[i]; //获取每一栋建筑
        var buildingHeight = currentBuilding.properties.height;
        var edgeNumber = currentBuilding.coordinates[0].length; //棱的个数

        turf.coordEach(buildings[i], function (currentCoord, coordIndex) {
            if (coordIndex == edgeNumber - 1) {
                nextCoord = currentBuilding.coordinates[0][0];
            } else {
                nextCoord = currentBuilding.coordinates[0][coordIndex + 1];
            }
            //顺时针旋转，北边是0°
            var currentAngle = turf.bearingToAzimuth(turf.bearing(observerCenter, currentCoord));//从前往后的向量,
            var nextAngle = turf.bearingToAzimuth(turf.bearing(observerCenter, nextCoord));

            //var point1 = turf.toMercator(Array.from(currentCoord));
            //var point2 = turf.toMercator(Array.from(nextCoord));
            //var point3 = turf.toMercator(Array.from(nextCoord));
            //var point4 = turf.toMercator(Array.from(currentCoord));

            //获取建筑物面的四个顶点
            var point1 = (Array.from(currentCoord)); var point2 = (Array.from(nextCoord));
            var point3 = (Array.from(nextCoord)); var point4 = (Array.from(currentCoord));
            point1.push(0);
            point2.push(0);
            point3.push(buildingHeight);
            point4.push(buildingHeight);
            shape = [point1, point2, point3, point4, point1];

            //如果立面不在同一直线上：将建筑物面的四个顶点加入到集合中
            if (currentAngle != nextAngle) {
                if (Math.abs(currentAngle - nextAngle) > 180) {
                    var shape = new buildingsShapes(shape, buildingHeight, [Math.max(currentAngle, nextAngle),
                    Math.min(nextAngle, currentAngle) + 360], [], []);
                    //console.log([Math.max(currentAngle, nextAngle),Math.min(nextAngle, currentAngle)+360])
                } else {
                    var shape = new buildingsShapes(shape, buildingHeight, [Math.min(currentAngle, nextAngle),
                    Math.max(nextAngle, currentAngle)], [], []);
                    //console.log([Math.min(currentAngle, nextAngle),Math.max(nextAngle, currentAngle)])
                }
                shapesGroup.push(shape);
            }
        });
    }
    //console.log(shapesGroup);
    return shapesGroup;
}

function isNaN(vald) {
    if (vald !== vald) {
        return true;
    }
    return false;
}

//计算直线
function calLine(p1, p2) {
    //https://www.cnblogs.com/pluslius/p/13800167.html
    //a = y2-y1, b = x1-x2, c = x2y1-x1y2
    var A = p2[1] - p1[1];
    var B = p1[0] - p2[0];
    var C = p2[0] * p1[1] - p1[0] * p2[1];
    //console.log(p1, p2);
    return [A, B, C];
}

//计算二维直线的交点
function calCross2DLine(l1, l2) {
    //x = (c2 * b1 - c1 * b2) / (a1 * b2 - a2 * b1)
    var x = (l2[2] * l1[1] - l1[2] * l2[1]) / (l1[0] * l2[1] - l2[0] * l1[1]);
    //y = (c1 * a2 - c2 * a1) / (a1 * b2 - a2 * b1)
    var y = (l1[2] * l2[0] - l2[2] * l1[0]) / (l1[0] * l2[1] - l2[0] * l1[1]);

    //console.log(x,y);
    var ll1 = l1[0] * x + l1[1] * y + l1[2];
    var ll2 = l2[0] * x + l2[1] * y + l2[2];

    if (isNaN(x) || isNaN(y)) //||ll1>0.1||ll2>0.1
    {
        //console.log(l1,l2);
        console.log(ll1, ll2, [x, y]);
        return NaN;
    }
    return [x, y];
}
//console.log(calCross2DLine([0,0,1], l2))

//计算单个面的遮挡面的遮挡面积
function calSingleShelterShape(currentShape, shelterShape, centerPoint) {
    /*
    currentShape：需要计算可视面积的面
    shelterShape：遮挡此面的面
    */

    //console.log(currentShape, shelterShape);
    var shape = [];
    var firstPoint;
    var minHeight = 1000;
    var nearCenter = false;

    for (var i = 0; i < 4; i++) {

        //获取延长线
        if (i < 2) {
            var l1 = calLine(shelterShape.shape[i], centerPoint); //可视中心与遮挡面构成的直线
            //console.log(centerPoint);
            var l2 = calLine(currentShape.shape[0], currentShape.shape[1]); //被遮挡面所构成的面
            //console.log(centerPoint,currentShape.shape[0], currentShape.shape[1]);
            //if (i == 0) {
            //如果两个射线角度小于90度，则相交点为这个点

            if (isNaN(calCross2DLine(l1, l2))) {
                //return [cenrerPoint, minHeight];
                var crossPoint = centerPoint;
                nearCenter = true;
            } else {
                var crossPoint = calCross2DLine(l1, l2);
                crossPoint.push(0);
                
                var v1 = [shelterShape.shape[i][0] -centerPoint[0], shelterShape.shape[i][1] - centerPoint[1]];
                var v2 = [crossPoint[0] - centerPoint[0], crossPoint[1] - centerPoint[1]];
                var angle = lineDotMultiply(v1, v2, option = "vector");
                if(angle<0){
                    
                    a = turf.distance(crossPoint,shelterShape.shape[i]);
                    b = turf.distance(crossPoint,shelterShape.shape[1-i]);
                    if (a > b){
                        crossPoint = Array.from(currentShape.shape[1-i]);
                    }else{
                        crossPoint = Array.from(currentShape.shape[i]);
                    }
                    //console.log(a,b,shelterShape.shape[i],crossPoint);
                    //有可能是crossPoint = Array.from(currentShape.shape[1-i]);
                    crossPoint.splice(2,1);
                    console.log(crossPoint);
                    crossPoint.push(0);
                }
                /**/
            }

            //console.log(calCross2DLine(l1, l2));
            

            shape.push(crossPoint);
            if (i == 0) {
                firstPoint = JSON.parse(JSON.stringify(crossPoint));
            }
        } else {
            //求高度 
            if (nearCenter) {
                var height = 10000;
            } else {
                var height = shelterShape.height / turf.distance(shelterShape.shape[3 - i], centerPoint) *
                    turf.distance(shape[3 - i], centerPoint);
            }

            var crossPoint = Array.from(shape[3 - i]);
            crossPoint[2] = height;
            shape.push(crossPoint);

            if (minHeight > height) {
                minHeight = height;
            }
        }
    }
    shape.push(firstPoint);
    return [shape, minHeight];
}

//计算遮挡建筑物面的面积
function calShelterBuilding(shapesGroup, observedArea) {
    console.log(shapesGroup);
    var point = (observedArea.observerCenter);//turf.toMercator
    point.push(0);

    outer:
    for (var i = 0; i < shapesGroup.length; i++) {
        //var angle = shapesGroup[i].angle;
        var shapes = [];

        //如果有遮挡面
        for (var j = 0; j < shapesGroup.length; j++) {
            if (i == j) continue;
            if (shapesGroup[i].isShelter(shapesGroup[j], point)) {
                //var currentCoord12 = shapesGroup[i].shape[0], currentCoord22 = shapesGroup[i].shape[1];

                //计算遮挡面投影到被遮挡面上得shape
                var [shape, minHeight] = calSingleShelterShape(shapesGroup[i], shapesGroup[j], point); //输入第一个需要计算阴影的面，第二是当前的投影面
                if (isNaN(shape)) {
                    continue;
                }
                shapes.push(shape);

                //加个判断条件，不符合条件的全部剔除,如果完全被影子笼罩则剔除
            
                if (shapesGroup[j].angle[0] <= shapesGroup[i].angle[0] && shapesGroup[j].angle[1] >= shapesGroup[i].angle[1]
                ) {//&& shapesGroup[i].height < minHeight
                    shapesGroup.splice(i, 1);
                    i--;
                    continue outer;
                }    /**/
            }
        }
        //console.log(shapes);
        shapesGroup[i].shelterShape = shapes;
    }
}

//向量叉乘计算：输入两个线段，线段的方向决定向量的方向
function lineCrossMultiply(l1, l2, option = "lineSegment") {
    if (option == "lineSegment") {
        //根据两直线计算：前-后
        var vectorL1 = [l1[0][0] - l1[1][0], l1[0][1] - l1[1][1], l1[0][2] - l1[1][2]];
        var vectorL2 = [l2[0][0] - l2[1][0], l2[0][1] - l2[1][1], l2[0][2] - l2[1][2]];
        // console.log(vectorL1, vectorL2);
        //计算叉乘
        var A = vectorL1[1] * vectorL2[2] - vectorL2[1] * vectorL1[2];
        var B = vectorL2[0] * vectorL1[2] - vectorL1[0] * vectorL2[2];
        var C = vectorL1[0] * vectorL2[1] - vectorL2[0] * vectorL1[1];

    } else if (option == "vector") {
        var A = l1[1] * l2[2] - l2[1] * l1[2];
        var B = l2[0] * l1[2] - l1[0] * l2[2];
        var C = l1[0] * l2[1] - l2[0] * l1[1];
    }
    return [A, B, C, A + B + C];
}

function lineDotMultiply(l1, l2, option = "lineSegment") {
    if (option == "lineSegment") {
        var vectorL1 = [l1[0][0] - l1[1][0], l1[0][1] - l1[1][1], l1[0][2] - l1[1][2]];
        var vectorL2 = [l2[0][0] - l2[1][0], l2[0][1] - l2[1][1], l2[0][2] - l2[1][2]];

        return vectorL1[0] * vectorL2[0] + vectorL1[1] * vectorL2[1] + vectorL1[2] * vectorL2[2];

    } else if (option == "vector") {
        return l1[0] * l2[0] + l1[1] * l2[1];
    }
}

function lineDotMultiply2(l1, l2, option = "lineSegment") {
    if (option == "lineSegment") {
        var vectorL1 = [l1[0][0] - l1[1][0], l1[0][1] - l1[1][1], l1[0][2] - l1[1][2]];
        var vectorL2 = [l2[0][0] - l2[1][0], l2[0][1] - l2[1][1], l2[0][2] - l2[1][2]];

        return vectorL1[0] * vectorL2[0] + vectorL1[1] * vectorL2[1] + vectorL1[2] * vectorL2[2];

    } else if (option == "vector") {
        return l1[0] * l2[0] + l1[1] * l2[1] + l1[2] * l2[2];
    }
}

//从总坐标转换为以建筑物为底面的坐标
function coordSystemConversion3To2(origin, currentShape) {
    //console.log("currentShape",currentShape);
    /*参数：原点，需要转换的面*/

    var shape2D = [];
    for (var i = 0; i < currentShape.length - 1; i++) {
        var currentPoint = currentShape[i];
        var distance2D = Math.sqrt(((currentPoint[0] - origin[0]) * (currentPoint[0] - origin[0]) +
            (currentPoint[1] - origin[1]) * (currentPoint[1] - origin[1])), 2); //求与原点的距离

        //计算两向量点乘大于零说明在同一个方向上
        if (lineDotMultiply2([currentPoint, origin], [currentShape[1], origin]) < 0) {
            distance2D *= -1;
        }
        shape2D.push([distance2D, currentPoint[2]]);
    }
    //console.log("shape2D",shape2D);
    shape2D.push(shape2D[0]);
    return shape2D;
}

//从以建筑物为底面的坐标转换为总坐标
function coordSystemConversion2To3(origin, coordSystem, shape2D) {


    coord = turf.getCoords(shape2D);
    //console.log(coord);
    coordSet = [];
    //for(var j = 0; j < coord.length;j++){
    if (coord.length == 1) {
        var coordShape = [];
        for (var i = 0; i < coord[0].length; i++) {
            //console.log(coord[j]);
            var currentShape = coord[0][i];
            var x = origin[0] + coordSystem[0] * (currentShape[0] / coordSystem[2]); //lon--x
            var y = origin[1] + coordSystem[1] * (currentShape[0] / coordSystem[2]); //lat--y
            coordWgs = turf.getCoords([x, y]);//turf.toWgs84

            coordWgs.push(currentShape[1]);
            coordShape.push(coordWgs);
        }
        coordSet.push(coordShape);
        //console.log(coordSet);
        return coordSet;
    } else {
        for (var j = 0; j < coord.length; j++) { //对每一个shape遍历
            var coordShape = [];
            for (var i = 0; i < coord[j][0].length; i++) {
                //console.log(coord[j][0]);
                var currentShape = coord[j][0][i];

                var x = origin[0] + coordSystem[0] * (currentShape[0] / coordSystem[2]); //lon--x
                var y = origin[1] + coordSystem[1] * (currentShape[0] / coordSystem[2]); //lat--y
                coordWgs = turf.getCoords(([x, y]));//turf.toWgs84
                //console.log(x,y,coordWgs);
                coordWgs.push(currentShape[1]);
                coordShape.push(coordWgs);
            }
            coordSet.push(coordShape);
        }
        return coordSet;
    }
}

//turf方位角规则：纬度方向上的北方是0度，顺时针经度向东为90度
function calShapesShadow(shapesGroup) {

    console.log("未投影区域计算的", shapesGroup);
    visualShapes = [];
    outer:
    for (var i = 0; i < shapesGroup.length; i++) {
        //计算每个shelter对应
        var currentShape = shapesGroup[i].shape;
        //console.log(currentShape);

        var shadow = [];
        var currentShadow;
        shadow.push(Array.from(currentShape)); //对于每一个面都有一个非遮挡区域，初始值为shape

        //把每个面转换坐标系console.log(origin)
        var coordSystemOrigin = currentShape[0]; //原点
        //var coordSystemDirection = turf.bearingToAzimuth(turf.bearing(currentShape[0], currentShape[1]));    //坐标轴方向
        var firstCoord = currentShape[0];
        var secondCoord = currentShape[1];
        var length = Math.sqrt(((secondCoord[0] - firstCoord[0]) * (secondCoord[0] - firstCoord[0]) +
            (secondCoord[1] - firstCoord[1]) * (secondCoord[1] - firstCoord[1])), 2);

        var coordSystem = [secondCoord[0] - firstCoord[0], secondCoord[1] - firstCoord[1], length]; //坐标轴系统参数
        //console.log(coordSystemOrigin,coordSystem);

        var currentShape2D = coordSystemConversion3To2(coordSystemOrigin, currentShape); //被遮挡面坐标系转换

        //console.log(currentShape2D);
        currentShape2D = turf.polygon([currentShape2D]);

        //console.log(shapesGroup[i].shelterShape.length);
        //对每一个阴影进行遍历
        for (var j = 0; j < shapesGroup[i].shelterShape.length; j++) {
            //console.log(shapesGroup[i]);
            currentShadow = shapesGroup[i].shelterShape[j];
            //console.log("aaaa",coordSystemOrigin,currentShadow,currentShape);
            var currentShadow2D = coordSystemConversion3To2(coordSystemOrigin, currentShadow);
            //console.log("currentShadow2D",currentShadow,currentShadow2D);
            currentShadow2D = turf.polygon([currentShadow2D]);

            //利用求差方法更新非遮挡面
            currentShape2D = turf.difference(currentShape2D, currentShadow2D);
            //console.log("currentShadow2D,currentShape2D", currentShadow2D, currentShape2D);

            if (currentShape2D == null) { //说明已经筛到最后没有了
                shapesGroup.splice(i, 1);
                i--;
                continue outer;
            }
        }
        visualShapes.push(coordSystemConversion2To3(coordSystemOrigin, coordSystem, currentShape2D));
    }
    //console.log(visualShapes);
    return visualShapes;
}

function secondaryScreening(visualShapes) {
    outer: for (var k = 0; k < visualShapes.length; k++) {
        var currentBuilding = visualShapes[k];
        //console.log(currentBuilding,currentBuilding.length);

        for (var i = 0; i < currentBuilding.length; i++) {
            var currentShape = currentBuilding[i];
            //console.log(currentShape);

            firstPoint = currentShape[0];
            //console.log(firstPoint,currentShape.length);

            var maxLon = 0,
                minLon = firstPoint[0];
            var maxLat = 0,
                minLat = firstPoint[1];
            var maxHeight = 0,
                minHeight = firstPoint[2];
            for (var j = 0; j < currentShape.length; j++) {
                var currentPoint = currentShape[j];
                if (maxLon < currentPoint[0]) maxLon = currentPoint[0];
                if (minLon > currentPoint[0]) minLon = currentPoint[0];
                if (maxLat < currentPoint[1]) maxLat = currentPoint[1];
                if (minLat > currentPoint[1]) minLat = currentPoint[1];
                if (maxHeight < currentPoint[2]) maxHeight = currentPoint[2];
                if (minHeight > currentPoint[2]) minHeight = currentPoint[2];
            }
            //console.log(maxLon,minLon,maxLat,minLat,maxHeight,minHeight);
            /*
            if (maxLon - minLon < 0.000001 || maxLat - minLat < 0.000001 || maxHeight - minHeight < 0.01) {
                currentBuilding.splice(i, 1);
                i--;
                if (currentBuilding.length == 0) {
                    visualShapes.splice(k, 1);
                    k--;
                    continue outer;
                }
            }*/
        }
    }
}