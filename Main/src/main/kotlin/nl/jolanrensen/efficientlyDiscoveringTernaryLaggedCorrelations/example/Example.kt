package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.example

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.WindowSkipping
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.WindowSkipping.numberOfWindowPlacementsSaved
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.WindowSkipping.removeStraightLineContinuity
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.randomF64FlatArray
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.simplification.RamerDouglasPeucker
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3ts
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3tsWithSkipping
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3tsWithSkipping.ContinuitySkippingType
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.toList
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.F64FlatArray

/**
 * This file is provided to give a couple of examples for how to use [StampPearson], [StampPearson3ts], and [StampPearson3tsWithSkipping].
 */


fun main() {

    stamp3tsWithSkippingExample()


}

/** Gives an example for STAMP-Pearson-3TS with 3 time series. */
fun stampPearson3tsExample() {

    // let's generate three random time series for this sample
    val tsA: F64FlatArray = randomF64FlatArray(100)
    val tsB: F64FlatArray = randomF64FlatArray(100)
    val tsC: F64FlatArray = randomF64FlatArray(100)

    // we need an instance of [StampPearson3ts]
    val stampPearson3ts = StampPearson3ts(maxArraySize = 100)

    val windowSize = 10

    // lag bound limits the distance between the windows
    val lagBound = 10

    // say the sliding means and standard deviations of one of the time series are precalculated somewhere,
    // then we can provide those to the function as well, to save calculations. If not provided, the function will simply
    // calculate them itself.
    val (timeSeriesBSlidingMeans: F64FlatArray, timeSeriesBSlidingStds: F64FlatArray) = StampPearson.computeSlidingMeanStd(
        tsB,
        windowSize)

    // if we want to debug and get all the correlation calculations, we can provide a 3D array to the function
    // this does make the calculation slower
    val resultArray = F64Array(
        tsC.length - windowSize + 1,
        tsB.length - windowSize + 1,
        tsA.length - windowSize + 1,
    )

    // let's call it!
    val bestCorrelation: Double = stampPearson3ts.stampPearson3ts(
        timeSeriesA = tsA,
        timeSeriesB = tsB,
        timeSeriesC = tsC,
        windowSize = windowSize,
        reducer = MAX,
        lagBound = lagBound,
        timeSeriesBSlidingMeans = timeSeriesBSlidingMeans,
        timeSeriesBSlidingStds = timeSeriesBSlidingStds,
        resultCubeCallback = resultArray,
    ) { aggregation: StampPearson3ts.AggregationMethod, aIndex: Int, bIndex: Int, cIndex: Int ->
        // we can provide a callback as last argument, which will give the indices and aggregation method of the best window positions found

        println("""
                |Best index in A is found at: $aIndex
                |Best index in B is found at: $bIndex
                |Best index in C is found at: $cIndex
                |Aggregation method used is: ${aggregation.textual}
            """.trimMargin())
    }

    println("best correlation found is $bestCorrelation")

    println(
        """ |all other correlation results are: 
            |${resultArray}""".trimMargin()
    )
}

/** Gives an example for STAMP-Pearson-3TS with 3 time series, simplification, and skipping. */
fun stamp3tsWithSkippingExample() {

    // let's generate three random time series for this sample
    val timeSeries: List<F64FlatArray> = buildList {
        repeat(3) {
            this += randomF64FlatArray(size = 100)
        }
    }

    val windowSize = 10

    // let's apply simplification to the time series using RamerDouglasPeucker
    val epsilon = 0.01
    val rdp = RamerDouglasPeucker(epsilon)
    timeSeries.forEach { rdp.filterKeepingLengthInPlace(it) }

    // Next, let's apply straight line skipping to a copy of the time series, while also keeping the original
    // We must also sort them such that tsA has the fewest NaN values and tsC the most
    val (tsA, tsB, tsC) = timeSeries
        .map {
            object {
                val original: F64FlatArray = it
                val withSkipping: F64FlatArray = it.copy().also { ts ->
                    removeStraightLineContinuity(ts, windowSize)
                }
            }
        }
        .sortedBy { numberOfWindowPlacementsSaved(it.withSkipping, windowSize) }


    // we need an instance of [StampPearson3tsWithSkipping]
    val stampPearson3tsWithSkipping = StampPearson3tsWithSkipping(
        maxArraySize = 100,
        continuitySkipping = ContinuitySkippingType.STRAIGHT_LINE,
    )


    // lag bound limits the distance between the windows
    val lagBound = 10

    // if we want to debug and get all the correlation calculations, we can provide a 3D array to the function
    // this does make the calculation slower
    val resultArray = F64Array(
        tsC.original.length - windowSize + 1,
        tsB.original.length - windowSize + 1,
        tsA.original.length - windowSize + 1,
    )

    // let's call it!
    val bestCorrelation: Double = stampPearson3tsWithSkipping.stampPearson3ts(
        timeSeriesA = tsA.withSkipping,
        timeSeriesB = tsB.withSkipping,
        timeSeriesC = tsC.withSkipping,
        windowSize = windowSize,
        reducer = MAX,

        timeSeriesAWithoutNaN = tsA.original,
        timeSeriesBWithoutNaN = tsB.original,
        timeSeriesCWithoutNaN = tsC.original,

        lagBound = lagBound,
        resultCubeCallback = resultArray,
    ) { aggregation: StampPearson3tsWithSkipping.AggregationMethod, aIndex: Int, bIndex: Int, cIndex: Int ->
        // we can provide a callback as last argument, which will give the indices and aggregation method of the best window positions found

        println("""
                |Best index in A is found at: $aIndex
                |Best index in B is found at: $bIndex
                |Best index in C is found at: $cIndex
                |Aggregation method used is: ${aggregation.textual}
            """.trimMargin())
    }

    println("best correlation found is $bestCorrelation")

    println(
        """ |all other correlation results are: 
            |${resultArray}""".trimMargin()
    )
}

/** Gives example of how to use STAMP-Pearson with 2 time series using this codebase. */
fun stampPearsonExample() {

    // let's generate two random time series for this sample
    val tsA: F64FlatArray = randomF64FlatArray(100)
    val tsB: F64FlatArray = randomF64FlatArray(100)

    // first we need an instance of [StampPearson]
    val stampPearson = StampPearson(maxArraySize = 100)

    val windowSize = 10

    // if we want to debug and get all the correlation calculations, we can provide a 2D array to the function
    // this does make the calculation slower
    val resultArray = F64Array(
        tsB.length - windowSize + 1,
        tsA.length - windowSize + 1,
    )

    val bestCorrelation: Double = stampPearson.stampPearson(
        timeSeriesA = tsA.copy(), // it's good practice providing a copy to the array, since they are likely overwritten by the algorithm
        timeSeriesB = tsB.copy(),
        windowSize = windowSize,
        reducer = MAX,
        resultPlaneCallback = resultArray
    ) { aIndex: Int, bIndex: Int ->
        // we can provide a callback as last argument, which will give the indices of the best window positions found

        println("""
                |Best index in A is found at: $aIndex
                |Best index in B is found at: $bIndex
            """.trimMargin())
    }

    println("best correlation found is $bestCorrelation")

    println(
        """ |all other correlation results are: 
            |${
            resultArray
                .toGenericArray()
                .map { (it as DoubleArray).toList() }
                .joinToString("\n")
        }""".trimMargin()
    )
}

/** Gives example of how to create a Matrix Profile using this codebase. */
fun matrixProfileExample() {

    // let's generate a random time series for the sample
    val timeSeries: F64FlatArray = randomF64FlatArray(100)

    // first we need an instance of [StampPearson]
    val stampPearson = StampPearson(maxArraySize = 100)

    val windowSize = 10

    // we'll call it using MAX correlation (so min Euclidean distance)
    val (matrixProfilePearson: F64FlatArray, matrixProfileIndices: IntArray) = stampPearson.matrixProfilePearson(
        timeSeries = timeSeries,
        windowSize = windowSize,
        reducer = MAX,
    )

    // if we want the Euclidean distance, we can simply convert the [matrixProfilePearson]
    val matrixProfile: F64FlatArray = StampPearson.pearsonCorrelationsToEuclidean(
        series = matrixProfilePearson,
        windowSize = windowSize,
    )

    // and we're done!
    println("""
        Using time series: ${timeSeries.toList()}
        We get the Matrix Profile: ${matrixProfile.toList()}
        With indices: ${matrixProfileIndices.toList()}
    """.trimIndent())
}