function angleToRadian(angle) {
    return angle * Math.PI / 180;
};

function radianToAngle(radian) {
    return radian * 180 / Math.PI;
};
//输入一段代码试试

//计算可视区域形状
function calVisualArea(brandCenterPoint, z, positionAngle, xResolution = 0.01, isAngle = true,
    eyeResolution = 3) {

    var x = brandCenterPoint[0],
        y = brandCenterPoint[1];
    //console.log("xResolution = ", z);
    if (isAngle == true) {
        eyeResolution = (eyeResolution / 60) / 60;
        eyeResolution = (eyeResolution * Math.PI) / 180; //人眼分辨率，弧度
    }

    var D = xResolution / eyeResolution;
    //半径
    var visualR = D / 2;
    var visualGroundR = Math.sqrt((Math.pow(D, 2)) / 4 - (Math.pow(z, 2))); //地面上的可视化半径
    //console.log("visualGroundR = ", Math.pow(z, 2));

    //坐标转换
    var visualHeight = z;
    visualR /= 1000;
    visualGroundR /= 1000;
    visualHeight /= 1000;
    brandCenterPoint = turf.toWgs84(brandCenterPoint);

    var visualCenter = turf.rhumbDestination(brandCenterPoint, visualGroundR,
        positionAngle, {
            units: 'kilometers'
        });
    //console.log("destination", visualCenter);


    var visualArea = {
            'brandCenterPoint': brandCenterPoint,
            'visualR': visualR,
            'visualGroundR': visualGroundR,
            'visualCenter': visualCenter,
            'visualHeight': visualHeight
        }
        //console.log("visualArea = ", visualArea);
    return visualArea;
}

//获取地面的可视坐标点
function getCirclePosition(visualArea) {

    //console.log(visualArea);
    var options = {
        steps: 360,
        units: 'kilometers'
    };

    //方法一：生成圆弧线转面
    //var circleLine = turf.lineArc(visualArea.visualCenter, visualArea.visualGroundR, 0, 360, options); //计算圆所有的点
    //console.log("circleLine = ",circleLine);
    //var circlePoly = turf.polygonize(circleLine); //线段转多边形
    //console.log(circleLine,"circlePoly = ",circlePoly);

    //生成圆弧面sector 
    var circlePoly = turf.circle(visualArea.visualCenter, visualArea.visualGroundR, options);
    //console.log("circlePoly = ", circlePoly);

    return circlePoly;
}

//获取可视建筑
function getVisualBuilding(circlePoly, buildings) {
    var total = 0;
    var visualBuilding = [];
    var idSelected = [];
    var exceptBuildingsPoly = circlePoly;

    //统计有多少个多边形相交
    turf.geomEach(buildings, function(currentGeometry, featureIndex, featureProperties) {
        if (turf.booleanOverlap(circlePoly, currentGeometry) || turf.booleanWithin(currentGeometry, circlePoly)) { //如果包含或者相交

            total++; //统计个数
            currentGeometry.properties = featureProperties; //特性移植
            //console.log(currentGeometry);

            //将面积连接

            var buildingShape = turf.difference(circlePoly, currentGeometry);
            try {
                exceptBuildingsPoly = turf.intersect(exceptBuildingsPoly, buildingShape);
            } catch (e) {
                //console.log(e);
            }

            //console.log(buildingsUnion);

            visualBuilding.push(currentGeometry); //记录可视范围内的建筑
            idSelected.push(currentGeometry.properties.id); //记录可视范围内的建筑编号
            //console.log(visualBuilding);
            //var line = turf.polygonToLine(currentGeometry); //

        };
    });


    var multiBuildPoly = turf.multiPolygon(visualBuilding); //多个多边形
    //console.log("idSelected = ", idSelected);
    //console.log("multiPoly = ", multiBuildPoly);
    //console.log("total = ", total);
    //console.log(buildingsUnion);

    var visualBuildings = {
        "multiBuildPoly": multiBuildPoly, //事业内的建筑
        "idSelected": idSelected, //视野内的建筑ID
        'exceptBuildingsPoly': exceptBuildingsPoly
    };

    return visualBuildings;
}

//计算阴影
function calShadow(circlePoly, visualBuildings, visualArea) {

    //对边界遍历获取阴影部分
    var adHeight = visualArea.visualHeight * 1000; //可视区域高度
    var adCoord = visualArea.brandCenterPoint; //广告牌的位置
    var visualCoord = visualArea.visualCenter.geometry.coordinates; //可视范围的中心点
    var visualGroundR = visualArea.visualGroundR; //可视区半径

    //console.log("visualArea", adCoord);
    var buildingHeight;
    var multiBuildPoly = visualBuildings.multiBuildPoly; //
    var buildNumber = multiBuildPoly.geometry.coordinates.length; //建筑物个数

    //console.log("multiBuildPoly = ", multiBuildPoly)

    //外包矩形

    var bbox = turf.bbox(circlePoly); //
    var circleR = Math.max(bbox[2] - bbox[0], bbox[3] - bbox[1]);

    //circleR *= 100;
    //console.log("circleR = ", circleR);

    var union = visualBuildings.exceptBuildingsPoly;

    var buildingShadow = [];

    for (var i = 0; i < buildNumber; i++) { //对每一栋建筑进行遍历buildNumber

        var currentBuilding = multiBuildPoly.geometry.coordinates[i]; //获取每一栋建筑
        //判断高度是否比广告牌高
        //console.log(currentBuilding);
        var buildingHeight = currentBuilding.properties.height;
        //console.log(buildingHeight);
        var edgeNumber = currentBuilding.coordinates[0].length; //棱的个数

        if (adHeight > (buildingHeight + 3)) { //如果广告牌高度更高

            turf.coordEach(currentBuilding, function(currentCoord, coordIndex) { //每条边
                //获取对应的点
                //console.log(currentCoord, coordIndex);//当前点
                let nextCoord
                    //给定下一个点的编号
                if (coordIndex == edgeNumber - 1) {
                    nextCoord = currentBuilding.coordinates[0][0];
                    //console.log(nextCoord);
                } else {
                    nextCoord = currentBuilding.coordinates[0][coordIndex + 1];
                    //console.log(nextCoord);
                }

                //先算第一个点的位置
                const currentX = -buildingHeight / (adHeight - buildingHeight) * (adCoord[0] - currentCoord[0]) + currentCoord[0];
                const currentY = -buildingHeight / (adHeight - buildingHeight) * (adCoord[1] - currentCoord[1]) + currentCoord[1];
                const nextX = -buildingHeight / (adHeight - buildingHeight) * (adCoord[0] - nextCoord[0]) + nextCoord[0];
                const nextY = -buildingHeight / (adHeight - buildingHeight) * (adCoord[1] - nextCoord[1]) + nextCoord[1];

                var polygon1 = turf.polygon([
                    [
                        [currentX, currentY],
                        [nextX, nextY], nextCoord, currentCoord, [currentX, currentY]
                    ]
                ], {
                    name: 'poly1'
                });
                //var polygon = turf.intersect(polygon1, circlePoly);
                //console.log(polygon2,polygon1);

                if (polygon1 != null) {

                    var polygonDiff = turf.difference(circlePoly, polygon1); //计算差异
                    if (polygonDiff != null && union != null) {
                        try {
                            union = turf.intersect(polygonDiff, union);
                        } catch (e) {
                            //console.log(e);
                        }
                    }
                    buildingShadow.push(polygon1);
                }
            });

        } else {
            turf.coordEach(currentBuilding, function(currentCoord, coordIndex) { //每条边
                //给定下一个点的编号
                let nextCoord
                if (coordIndex == edgeNumber - 1) {
                    nextCoord = currentBuilding.coordinates[0][0];
                    //console.log(nextCoord);
                } else {
                    nextCoord = currentBuilding.coordinates[0][coordIndex + 1];
                    //console.log(nextCoord);
                }

                var currentAngle = (turf.bearing(adCoord, currentCoord));//turf.bearingToAzimuth
                var nextAngle = (turf.bearing(adCoord, nextCoord));
                
                if (currentAngle != nextAngle) {

                    var arc = turf.lineArc(adCoord, visualGroundR * 2, turf.bearingToAzimuth(Math.min(currentAngle, nextAngle)),
                        turf.bearingToAzimuth(Math.max(nextAngle, currentAngle))).geometry.coordinates;
                        
                    arc.push(nextCoord);
                    arc.unshift(currentCoord);
                    arc.push(currentCoord);

                    var polygon = turf.polygon([arc]);
                    //var polygon = turf.intersect(polygon1, circlePoly);
                    //console.log(polygon,polygon1);
                    if (polygon != null) {

                        var polygonDiff = turf.difference(circlePoly, polygon); //计算差异
                        if (polygonDiff != null && union != null) {

                            //console.log(buildingHeight, polygonDiff, union);
                            try {
                                union = turf.intersect(polygonDiff, union);
                            } catch (e) {
                                //console.log(e);
                            }

                            //.catch((e) => {});
                        }

                        buildingShadow.push(polygon);
                    }

                }

            }); //coordEach
        } //else
    } //for building

    var buildingShadow = turf.featureCollection(buildingShadow);
    let visualPoly
    if (union != null) {
        visualPoly = turf.intersect(circlePoly, union);
    } else {
        visualPoly = null;
    }

    visualAreaShadow = {
        'buildingShadow': buildingShadow,
        'visualPoly': visualPoly

    }
    return (visualAreaShadow);
} //function