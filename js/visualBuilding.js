function angleToRadian(angle) {
    return angle * Math.PI / 180;
};

function radianToAngle(radian) {
    return radian * 180 / Math.PI;
};

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
    if(visualR>z){
        var visualGroundR = Math.sqrt((Math.pow(D, 2)) / 4 - (Math.pow(z, 2))); //地面上的可视化半径
    }else{
        var visualGroundR = 0;
    }
    
    //console.log("visualGroundR = ", Math.pow(z, 2));

    //坐标转换
    var visualHeight = z;
    visualR /= 1000;
    visualGroundR /= 1000;
    visualHeight /= 1000;
    brandCenterPoint = turf.toWgs84(brandCenterPoint);

    var visualCenter = turf.rhumbDestination(brandCenterPoint, visualGroundR,
        positionAngle, { units: 'kilometers' });
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
    var options = { steps: 360, units: 'kilometers' };

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

//获取可视范围内的建筑
function getVisualBuilding(circlePoly, buildings) {
    var total = 0;
    var visualBuilding = [];
    var idSelected = [];
    //var exceptBuildingsPoly = circlePoly;
    
    var start = new Date().getTime();
    //统计有多少个多边形相交
    turf.geomEach(buildings, function(currentGeometry,featureIndex,  featureProperties) {
        if (turf.booleanOverlap(circlePoly, currentGeometry) || turf.booleanWithin(currentGeometry, circlePoly)) { //如果包含或者相交
            var start = new Date().getTime();
            total++; //统计个数
            currentGeometry.properties = featureProperties; //特性移植
            //console.log(currentGeometry);

            visualBuilding.push(currentGeometry); //记录可视范围内的建筑
            idSelected.push(currentGeometry.properties.id); //记录可视范围内的建筑编号
            var end = new Date().getTime();
            console.log('cost is', `${end - start}ms`);
            //console.log(currentGeometry,featureProperties);
            //var line = turf.polygonToLine(currentGeometry); //
        };
    });
   


    var multiBuildPoly = turf.multiPolygon(visualBuilding); //多个多边形

    var visualBuildings = {
        "multiBuildPoly": multiBuildPoly, //事业内的建筑
        "idSelected": idSelected, //视野内的建筑ID
        //'exceptBuildingsPoly': exceptBuildingsPoly
    };
    
    return visualBuildings;
}

//获取可见建筑物
function calVisibleBuilding(circlePoly, visualBuildings, visualArea){
    
}



//计算阴影
function calBuildingsShadow(circlePoly, visualBuildings, visualArea) {

    //对边界遍历获取阴影部分
    var adHeight = visualArea.visualHeight * 1000; //可视区域高度
    var adCoord = visualArea.brandCenterPoint; //广告牌的位置
    //var visualCoord = visualArea.visualCenter.geometry.coordinates; //可视范围的中心点
    var visualGroundR = visualArea.visualGroundR; //可视区半径

    //console.log("visualArea", adCoord);
    var buildingHeight;
    var multiBuildPoly = visualBuildings.multiBuildPoly; //
    var buildNumber = multiBuildPoly.geometry.coordinates.length; //建筑物个数

    //console.log("multiBuildPoly = ", multiBuildPoly)

    var buildingShadow = [];

    for (var i = 0; i < buildNumber; i++) { //对每一栋建筑进行遍历buildNumber

        var currentBuilding = multiBuildPoly.geometry.coordinates[i]; //获取每一栋建筑
        //判断高度是否比广告牌高
        //console.log(currentBuilding);
        var buildingHeight = currentBuilding.properties.height;
        //console.log(buildingHeight);
        var edgeNumber = currentBuilding.coordinates[0].length; //棱的个数

        buildingShadow.push(turf.polygon(currentBuilding.coordinates))

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

                var polygon = turf.polygon([
                    [
                        [currentX, currentY],
                        [nextX, nextY], nextCoord, currentCoord, [currentX, currentY]
                    ]
                ], { name: 'poly1' });

                if (polygon != null) {
                    buildingShadow.push(polygon);
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
                }

                var currentAngle = turf.bearingToAzimuth(turf.bearing(adCoord, currentCoord)); //turf.bearingToAzimuth
                var nextAngle = turf.bearingToAzimuth(turf.bearing(adCoord, nextCoord));

                if (currentAngle != nextAngle) {
                    if (Math.abs(currentAngle - nextAngle) > 180) {
                        var arc = turf.lineArc(adCoord, visualGroundR * 2, Math.max(currentAngle, nextAngle),
                            Math.min(nextAngle, currentAngle)).geometry.coordinates;
                    } else {
                        var arc = turf.lineArc(adCoord, visualGroundR * 2, Math.min(currentAngle, nextAngle),
                            Math.max(nextAngle, currentAngle)).geometry.coordinates;
                    }

                    //console.log(nextCooord,currentCoord)
                    arc.push(nextCoord);
                    arc.unshift(currentCoord);
                    arc.push(currentCoord);

                    var polygon = turf.polygon([arc]);
                    //var polygon = turf.intersect(polygon1, circlePoly);
                    //console.log(polygon,polygon1);
                    if (polygon != null) {

                        buildingShadow.push(polygon);
                    }
                }

            }); //coordEach
        } //else
    } //for building

    var buildingShadow = turf.featureCollection(buildingShadow);

    const visualAreaShadow = {
        'buildingShadow': buildingShadow,
        //'visualPoly': visualPoly
    }
    return (visualAreaShadow);
} //function

//计算最后的可视区域
function calvisibleArea(circlePoly, buildingShadow) {
    var polygonDiff;
    var union = circlePoly;
    //统计建筑物地面形状
    turf.geomEach(buildingShadow, function(currentGeometry, featureIndex) {

        //console.log(currentGeometry, featureIndex);
        polygonDiff = turf.difference(circlePoly, currentGeometry); //计算差异
        if (polygonDiff != null && union != null) {
            //console.log(union);
            try {
                union = turf.intersect(polygonDiff, union);
            } catch (e) {
                //console.log(e);
            }
        }
    });
    //console.log(union);
    return union;
}