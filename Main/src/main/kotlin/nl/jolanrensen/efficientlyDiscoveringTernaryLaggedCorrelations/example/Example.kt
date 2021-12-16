package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.example

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.randomF64FlatArray
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.toList
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor.asF64Array

/**
 * This file is provided to give a couple of examples for how to use [StampPearson], [StampPearson3ts], and [StampPearson3tsWithSkipping].
 */


fun main() {

    // Let's generate some random time series first
    // Usually data is read from a database
    val tsA: F64FlatArray = randomF64FlatArray(100)
    val tsB: F64FlatArray = randomF64FlatArray(100)
    val tsC: F64FlatArray = randomF64FlatArray(100)

    // We use [F64FlatArrays] a lot, since they have memory efficient operations. [DoubleArray]s can also be converted
    val exampleTimeSeries: F64FlatArray = doubleArrayOf(1.0, 2.0, 3.0).asF64Array()


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
        timeSeriesA = tsA,
        timeSeriesB = tsB,
        windowSize = windowSize,
        reducer = MAX,
        resultPlaneCallback = resultArray,
        resultIndexCallback = { aIndex, bIndex ->
            println("""
                |Best index in A is found at: $aIndex
                |Best index in B is found at: $bIndex
            """.trimMargin())
        }
    )

    println("best correlation found is $bestCorrelation")

    println("all other correlation results are: $resultArray")
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