package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations

import nl.jolanrensen.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson
import org.jetbrains.bio.viktor.F64FlatArray
import org.nield.kotlinstatistics.simpleRegression
import java.io.Serializable
import java.util.*
import kotlin.math.abs

object Continuity : Serializable {

    class MotifInstance(val data: List<Int>) {
        val range = data.first()..data.last()
        override fun toString(): String = data.toString()
    }

    tailrec fun find(parent: IntArray, i: Int): Int = when {

        parent[i] == i || parent[i] == -1 ->
            i

        else ->
            find(parent, parent[i])
    }

    fun union(parent: IntArray, x: Int, y: Int) {
        require(x != y)
        val xSet = find(parent, x)
        val ySet = find(parent, y)

        parent[xSet] = ySet
    }

    fun removeRepeatingPatterns(
        timeSeries: F64FlatArray,
        windowSize: Int,
        matrixProfilePrecision: Double = 0.9,
        stampPearson: StampPearson,
    ) {
        require(timeSeries.length >= 2 && windowSize >= 2)

        val (matrixProfilePearson, indices) = stampPearson.matrixProfilePearson(
            timeSeries = timeSeries,
            windowSize = windowSize,
            reducer = MAX,
        )

//        val stage1 = timeSeries.copy()
//        println("matrix profile: ${matrixProfilePearson.toList()}")
//        println("indices mp: ${indices.toList()}")

        val motifInstances = mutableListOf<MotifInstance>()

        var currentMotifInstance: MutableList<Int> = mutableListOf()
        for (i in matrixProfilePearson.indices) {
            if (matrixProfilePearson[i] >= matrixProfilePrecision && (currentMotifInstance.isEmpty() || currentMotifInstance.last() == i - 1)) {
                currentMotifInstance.add(i)
            } else if (currentMotifInstance.isNotEmpty()) {
                motifInstances += MotifInstance(currentMotifInstance)
                currentMotifInstance = mutableListOf()
            }
        }

        // filter out the ones too short
        motifInstances.removeAll {
            when {
                0 in it.range || timeSeries.length - 1 in it.range -> it.data.size <= windowSize
                else -> it.data.size <= 2 * windowSize - 1
            }
        }
//        println(motifInstances)



        val clusterIndices = IntArray(motifInstances.size) { -1 }

        val matchesWithCount = mutableMapOf<Pair<Int, Int>, Int>()

//        data class DirectedEdge(val origin: Int, val destination: Int)

//        val edges = mutableSetOf<DirectedEdge>()

        for ((i, motifInstance) in motifInstances.withIndex()) {
            val pointingAt = motifInstance.data.map { indices[it] }

            data class Count(val j: Int, val count: Int)

            val topX = PriorityQueue<Count>(1) { a, b -> a.count.compareTo(b.count) }

            for ((j, otherInstance) in motifInstances.withIndex()) {
                if (i == j) continue // while they can point at themselves, this provides no information
                if (motifInstance.data.size != otherInstance.data.size) continue

                val count = pointingAt.count { it in otherInstance.range }

                if (count > 0) {
//                    edges += DirectedEdge(i, j)
                    topX += Count(j, count)
                }
            }

            topX.forEach { (j, _) ->
//                val key = minOf(i, j) to maxOf(i, j)
//                matchesWithCount[key] = (matchesWithCount[key] ?: 0) + 1
                union(parent = clusterIndices, x = i, y = j)
            }
        }

        // TODO here
//        for ((pair, count) in matchesWithCount) {
//            val (i, j) = pair
//            union(clusterIndices, i, j)
//            union(clusterIndices, j, i)
//        }


        val clusters = clusterIndices
            .filterNot { it == -1 }
            .toSet()
            .associateWith { mutableListOf<MotifInstance>() }

        for ((i, motifInstance) in motifInstances.withIndex())
            clusters[clusterIndices[i]]?.add(motifInstance)

        val motifs: List<List<MotifInstance>> = clusters.values.toList().onEach {
            // make sure the pieces where the largest pieces can be removed are given priority
            it.sortBy { it.data.first() == 0 || it.data.last() == timeSeries.length - 1 }
        }

        for (motif in motifs) {
            for ((i, motifInstance) in motif.withIndex()) {
                if (i == 0) continue

                // keep only first m - 1 and last m of slope unless at the end of the timeseries, then only keep first m
                if (motifInstance.data.last() == timeSeries.length - 1) { // TODO
                    // keep only first m - 1 of motif
                    motifInstance.data.takeLast(motifInstance.data.size - windowSize + 1).forEach {
//                        if (timeSeries[it].isNaN()) throw Exception("$it is being set to NaN multiple times")
                        timeSeries[it] = Double.NaN
                    }
                } else {
                    // keep only first m - 1 and last m - 1 of motif
                    motifInstance.data.slice(windowSize - 1..motifInstance.data.size - windowSize).forEach {
//                        if (timeSeries[it].isNaN()) throw Exception("$it is being set to NaN multiple times")
                        timeSeries[it] = Double.NaN
                    }
                }

                if (motifInstance.data.first() == 0) {
                    // keep only last m - 1 of motif
                    motifInstance.data.take(motifInstance.data.size - windowSize + 1).forEach {
//                        if (timeSeries[it].isNaN()) throw Exception("$it is being set to NaN multiple times")
                        timeSeries[it] = Double.NaN
                    }
                }

            }
        }

//        println("time series after: ${timeSeries.toList()}")
//        println("window size: $windowSize, removed ${timeSeries.count { it.isNaN() }}/${timeSeries.size} values. Saved ${
//            numberOfWindowPlacementsSaved(timeSeries,
//                windowSize)
//        }/${timeSeries.size - windowSize} window placements")

//        if (timeSeries.count { it.isNaN() } > 0)
//        quicklyPlotTimeSeries(stage1, matrixProfilePearson)
//        whileTrue()
    }


    @Suppress("ReplaceRangeToWithUntil")
    fun removeContinuity2(
        timeSeries: F64FlatArray,
        windowSize: Int,
        precision: Double = 0.01,
        rSquarePrecision: Double = 0.7,
    ) {
        require(timeSeries.length >= 2 && windowSize >= 2)

        var slopeStart: Int? = null
        var slopeEnd: Int? = null // inclusive
        var currentDelta: Double? = null

        val applyToTimeSeries = mutableMapOf<Int, Double>()

        timeSeries.mapWindowedWithIndex(windowSize) { start, query ->
            val end = start + windowSize - 1 // inclusive

            val regression = query
                .toDoubleArray()
                .mapIndexed { i, data -> start + i to data }
                .simpleRegression()

//            println("query: $query")
//            println("regression result: ${
//                query.toDoubleArray().indices.map { i -> regression.predict((start + i + 1).toDouble()) }
//            }")

            val slope = regression.slope.coerceNotNaN(replaceWith = Double.POSITIVE_INFINITY)
            val rSquare = regression.rSquare.coerceNotNaN(replaceWith = 0.0)


            fun startSlope() {
                slopeStart = start
                slopeEnd = end
                currentDelta = slope
            }

            fun endSlope() {
                // TODO fill in the last (or first) m with the values gained from the average slope

                // slope cannot continue
                val slopeSize = slopeEnd!! - slopeStart!!

                if (slopeSize > windowSize) {

                    if (slopeEnd != timeSeries.length - 1) {
                        // replace first m - 1 with average of slope TODO
                        for (j in slopeStart!!..(slopeStart!! + windowSize - 2)) {
                            applyToTimeSeries[j] = regression.predict(j.toDouble())
                        }

                        // keep only first m - 1 and last m of slope unless at the end of the time series
                        for (j in (slopeStart!! + windowSize - 1)..(slopeEnd!! - windowSize)) {
                            applyToTimeSeries[j] = Double.NaN
                        }

                        // replace last m of slope with queries of slope
                        for (j in (slopeEnd!! - windowSize + 1)..slopeEnd!!) {
                            applyToTimeSeries[j] = regression.predict(j.toDouble())
                        }
                    } else { // end of time series may be NaN too
                        for (j in (slopeStart!! + windowSize)..slopeEnd!!) {
                            applyToTimeSeries[j] = Double.NaN
                        }
                    }

                    // start of time series may be NaN too
                    if (slopeStart == 0) {
                        for (j in 0..(windowSize - 2)) { // the rest is done using next for loop
                            applyToTimeSeries[j] = Double.NaN
                        }
                    }


                }
            }

            when {

                // if current piece of a slope has right r²
                rSquare >= rSquarePrecision -> {

                    when {

                        // if at start
                        slopeStart == null || slopeEnd == null || currentDelta == null -> {
                            // start new slope
                            startSlope()
                        }

                        // try to continue previous slope
                        abs(currentDelta!! - slope) <= precision -> {
                            slopeEnd = end

                            // if at end of time series
                            if (end == timeSeries.length - 1) {
                                endSlope()
                            }
                        }

                        // else previous slope is too different, end and start new one
                        else -> {
                            endSlope()

                            // start new slope
                            startSlope()
                        }
                    }
                }

                // if r² is too large and slope is started
                slopeStart != null && slopeEnd != null && currentDelta != null -> {
                    endSlope()

                    // start new slope
                    startSlope()
                }

                // else no slope is started yet
                else -> Unit
            }
        }

        for ((i, double) in applyToTimeSeries) timeSeries[i] = double

        println("removed ${applyToTimeSeries.size} items")
    }

    /**
     * Removes continuity of slope larger than [windowSize] from the timeSeries by replacing it with [Double.NaN].
     *
     * If slope touches start, keep only m of end of slope
     * If slope touches end, keep only m of the start of the slope
     * Else keep m - 1 of start and m of end of slope
     *
     *
     *
     */
    @Suppress("ReplaceRangeToWithUntil")
    fun removeStraightLineContinuity(
        timeSeriesSrc: F64FlatArray,
        windowSize: Int,
        timeSeriesDest: F64FlatArray = timeSeriesSrc,
        precision: Double = 0.001,
    ) {
        require(timeSeriesSrc.length >= 2 && windowSize >= 2)

        var slopeStart = 0
        var slopeEnd = 1 // inclusive
        var slopeDelta = timeSeriesSrc[0] - timeSeriesSrc[1]
        var intermediatePrecision = precision

        for (i in 1 until timeSeriesSrc.length) {
            val delta = if (i == timeSeriesSrc.length - 1) Double.POSITIVE_INFINITY // for last round
            else timeSeriesSrc[i] - timeSeriesSrc[i + 1]

            if (abs(delta - slopeDelta) <= intermediatePrecision) {
                // slope continues
                slopeEnd = i + 1
                intermediatePrecision /= 2.0
            } else {
                // slope ends

                val slopeSize = slopeEnd - slopeStart

                if (slopeSize > windowSize) {

                    // keep only first m - 1 and last m of slope unless at the end of the timeseries
                    if (slopeEnd == timeSeriesSrc.length - 1) { // end of time series may be NaN too
                        for (j in (slopeStart + windowSize)..slopeEnd) {
                            if (timeSeriesDest[j].isNaN()) println("overwriting NaN")

                            timeSeriesDest[j] = Double.NaN
                        }
                    } else {
                        for (j in slopeStart..(slopeStart + windowSize - 2)) {
                            if (timeSeriesDest[j].isNaN()) println("overwriting NaN")

                            timeSeriesDest[j] = timeSeriesSrc[j]
                        }
                        for (j in (slopeStart + windowSize - 1)..(slopeEnd - windowSize)) {
                            if (timeSeriesDest[j].isNaN()) println("overwriting NaN")

                            timeSeriesDest[j] = Double.NaN
                        }
                        for (j in (slopeEnd - windowSize + 1)..slopeEnd) {
                            if (timeSeriesDest[j].isNaN()) println("overwriting NaN")

                            timeSeriesDest[j] = timeSeriesSrc[j]
                        }
                    }

                    // start of time series may be NaN too
                    if (slopeStart == 0) {
                        for (j in slopeStart..(slopeStart + windowSize - 2)) {
                            if (timeSeriesDest[j].isNaN()) println("overwriting NaN")

                            timeSeriesDest[j] = Double.NaN
                        }
                    }
                }

                // start new slope
                slopeStart = i
                slopeEnd = i + 1
                slopeDelta = delta
                intermediatePrecision = precision
            }
        }


        // TODO REMOVE
//        val noRemoved = timeSeriesDest.asIterable().count { it.isNaN() }
//        println("removed $noRemoved continuity items of ${timeSeriesDest.size}, precision: $precision, m: $windowSize")
    }

    fun numberOfWindowPlacementsSaved(timeSeriesWithNaNs: F64FlatArray, windowSize: Int): Int {
        var sum = 0

        var current = 0
        for (value in timeSeriesWithNaNs) {
            when {
                value.isNaN() -> current++
                current > 0 -> {
                    sum += windowSize + current - 1
                    current = 0
                }
            }
        }

        return sum
    }

    @Suppress("ReplaceRangeToWithUntil")
    @Deprecated("won't be possible for repeating patterns")
    fun restoreContinuity(timeSeries: F64FlatArray) {

        val n = timeSeries.length
        var noStartingNaNs = 0
        var noEndingNaNs = 0


        // if series starts with NaN
        if (timeSeries[0].isNaN()) {
            var i = 0
            while (timeSeries[i].isNaN()) i++

            noStartingNaNs = i
            // println("noStartingNaNs: $noStartingNaNs")

            while (timeSeries[i].isNotNaN() && i - noStartingNaNs + 1 < 2 && i < n - 1) i++

            val noValuesAfterNaNs = i - noStartingNaNs + 1
            // println("noValuesAfterNans: $noValuesAfterNaNs")

            require(noValuesAfterNaNs >= 2) { "Not enough values after starting NaNs to restore continuity." }
            for (j in noStartingNaNs - 1 downTo 0) {
                timeSeries[j] = timeSeries[j + 1] - (timeSeries[j + 2] - timeSeries[j + 1])
            }
        }

        // if series ends with NaN
        if (timeSeries[n - 1].isNaN()) {
            var i = n - 1
            while (timeSeries[i].isNaN()) i--

            noEndingNaNs = n - i - 1
            // println("noEndingNaNs: $noEndingNaNs")

            while (timeSeries[i].isNotNaN() && n - noEndingNaNs - i < 2 && i > 0) i--

            val noValuesBeforeNaNs = n - noEndingNaNs - i
            // println("noValuesBeforeNaNs: $noValuesBeforeNaNs")

            require(noValuesBeforeNaNs >= 2) { "Not enough values before ending NaNs to restore continuity." }

            for (j in n - noEndingNaNs..n - 1) {
                timeSeries[j] = timeSeries[j - 1] + (timeSeries[j - 1] - timeSeries[j - 2])
            }
        }


        // middle part
        for (i in n - noEndingNaNs - 1 downTo noStartingNaNs) {
            if (timeSeries[i].isNotNaN()) continue
            require(timeSeries[i + 1].isNotNaN() && timeSeries[i + 2].isNotNaN()) { "Not enough values in middle to restore continuity." }
            timeSeries[i] = timeSeries[i + 1] - (timeSeries[i + 2] - timeSeries[i + 1])
        }
    }
}