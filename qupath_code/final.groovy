import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.PathDetectionObject
import qupath.lib.geom.Point2
print "✅ 细胞检测完成"

// ✅ 2. 获取所有Annotation芯
def annotations = getAnnotationObjects()

// ✅ 3. 设置通道名（和你的实际数据对应）
def channelCDH1 = "Cell: Opal 480 max"   
def channelPanCK = "Cell: Opal 520 max"
def channelTIMP3 = "Cytoplasm: Opal 570 max"

// ✅ 4. 设置筛选阈值
def thresholdCDH1 = 100
def thresholdPanCK = 100
def thresholdTIMP3 = 100

// ✅ 5. 准备存储每组的距离列表
def groupDistances = [:] // key = group name ("sensitive" / "resistance"), value = list of distances

// ✅ 6. 遍历每个芯Annotation，按组别处理
annotations.each { annotation ->
    def groupName = annotation.getPathClass()?.toString()?.toLowerCase()
    if (groupName != "sensitive" && groupName != "resistance") {
        // 不是感兴趣的组别，跳过
        return
    }
    
    // 只找芯内部的细胞
    def cellsInAnnotation = getDetectionObjects().findAll { cell -> 
        cell.getROI() != null && annotation.getROI().contains(cell.getROI().getCentroidX(), cell.getROI().getCentroidY())
    }
    
    // 先把所有细胞打成 Non-Epi
    cellsInAnnotation.each { cell ->
        cell.setPathClass(getPathClass("Non-Epi"))
    }
    
    // 根据CDH1/PanCK表达重新打Tumor类
    cellsInAnnotation.each { cell ->
        def valueCDH1 = cell.getMeasurementList().get(channelCDH1)
        def valuePanCK = cell.getMeasurementList().get(channelPanCK)
        if ((valueCDH1 != null && valueCDH1 > thresholdCDH1) || (valuePanCK != null && valuePanCK > thresholdPanCK)) {
            cell.setPathClass(getPathClass("Tumor"))
        }
    }
    
    // 再筛选 Non-Epi中TIMP3+的基质细胞
    def timp3PositiveStroma = []
    cellsInAnnotation.each { cell ->
        if (cell.getPathClass() == getPathClass("Non-Epi")) {
            def valueTIMP3 = cell.getMeasurementList().get(channelTIMP3)
            if (valueTIMP3 != null && valueTIMP3 > thresholdTIMP3) {
                cell.setPathClass(getPathClass("Stroma_TIMP3+"))
                timp3PositiveStroma << cell
            }
        }
    }
    
    // 找到Tumor细胞
    def tumorCells = cellsInAnnotation.findAll { it.getPathClass() == getPathClass("Tumor") }
    
    if (tumorCells.isEmpty()) {
        print "⚠️ 组 ${groupName} 中没有Tumor细胞，跳过距离计算！"
        return
    }
    
    // 计算每个TIMP3+细胞到最近Tumor细胞的距离
    def distances = []
    for (cell in timp3PositiveStroma) {
        def roi = cell.getROI()
        def cellCenter = new Point2(roi.getCentroidX(), roi.getCentroidY())
        
        def minDistance = Double.MAX_VALUE
        for (tumor in tumorCells) {
            def tumorCenter = new Point2(tumor.getROI().getCentroidX(), tumor.getROI().getCentroidY())
            def distance = cellCenter.distance(tumorCenter)
            if (distance < minDistance) {
                minDistance = distance
            }
        }
        
        distances << minDistance
        // 也可以保存到细胞MeasurementList（可选）
        cell.getMeasurementList().putMeasurement("Distance_to_Nearest_Tumor", minDistance)
    }
    
    // 保存当前组的距离
    if (!groupDistances.containsKey(groupName)) {
        groupDistances[groupName] = []
    }
    groupDistances[groupName].addAll(distances)
}

// ✅ 7. 统计每组的均值、中位数、标准差、最小最大，并打印
groupDistances.each { group, distances ->
    if (distances.isEmpty()) {
        println "⚠️ 组 ${group} 没有找到TIMP3+基质细胞或没有计算到距离！"
        return
    }
    
    distances.sort()
    def n = distances.size()
    def mean = distances.sum() / n
    def sd = Math.sqrt(distances.collect { (it - mean) * (it - mean) }.sum() / (n - 1))
    def median = (n % 2 == 0) ? (distances[(n/2)-1] + distances[n/2])/2.0 : distances[(int)(n/2)]
    def min = distances.first()
    def max = distances.last()

    println "\n=== 组别: ${group.toUpperCase()} ==="
    println "细胞数目 (N)        = ${n}"
    println "平均距离 (Mean)     = ${String.format('%.2f', mean)} μm"
    println "标准差 (SD)         = ${String.format('%.2f', sd)} μm"
    println "中位数 (Median)     = ${String.format('%.2f', median)} μm"
    println "最小距离 (Min)      = ${String.format('%.2f', min)} μm"
    println "最大距离 (Max)      = ${String.format('%.2f', max)} μm"
}

fireHierarchyUpdate()

import qupath.lib.objects.PathDetectionObject
import qupath.lib.objects.classes.PathClass
import qupath.lib.scripting.QP
import java.io.FileWriter
import java.io.BufferedWriter

// ✅ 创建保存csv文件的路径（自己改成你想要的目录）
def outputFile = buildFilePath(PROJECT_BASE_DIR, 'exported_distances.csv')
def writer = new BufferedWriter(new FileWriter(outputFile))

// ✅ 写表头
writer.write("group,distance\n")

// ✅ 遍历每组，写每个细胞的 group 和 distance
groupDistances.each { group, distances ->
    distances.each { distance ->
        writer.write("${group},${distance}\n")
    }
}

writer.close()

print "✅ 已成功保存细胞group和distance到CSV：${outputFile}"

