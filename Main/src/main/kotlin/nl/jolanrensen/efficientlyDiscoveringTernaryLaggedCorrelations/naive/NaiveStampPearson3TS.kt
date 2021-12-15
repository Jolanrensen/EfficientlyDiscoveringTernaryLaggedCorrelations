package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.naive

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.AggregationWithReducerWithArgWithCompare
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.*
import org.jetbrains.bio.viktor.F64FlatArray
import java.io.Serializable
import kotlin.math.abs

object NaiveStampPearson3TS : Serializable {

    // lagbound
    fun pearson3dCube(
        timeSeriesA: F64FlatArray,
        timeSeriesB: F64FlatArray,
        timeSeriesC: F64FlatArray,
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,
        lagBound: Int? = null,
    ): Double {
        val m = windowSize
        val nC = timeSeriesC.length

        if (lagBound != null) {

            var best = reducer.initializer
            timeSeriesC.forEachWindowWithIndex(m) { queryCIndex, queryC ->
                val result = pearson3dPlane(
                    timeSeriesA = timeSeriesA,
                    timeSeriesB = timeSeriesB,
                    queryC = queryC,
                    aggReducer = reducer,
                    lagBound = lagBound,
                    queryCIndex = queryCIndex,
                )
                if (reducer.firstIsBetterThanSecond(first = result, second = best)) {
                    best = reducer(result, best)
                }
            }

            return best

        } else {
            var best = reducer.initializer

            timeSeriesC.forEachWindow(windowSize) { queryC ->
                val result = pearson3dPlane(
                    timeSeriesA = timeSeriesA,
                    timeSeriesB = timeSeriesB,
                    queryC = queryC,
                    aggReducer = reducer,
                )
                if (reducer.firstIsBetterThanSecond(first = result, second = best)) {
                    best = reducer(result, best)
                }
            }

            return best

        }
    }


    // 2 time series, 1 query, lagbound
    private fun pearson3dPlane(
        timeSeriesA: F64FlatArray,
        timeSeriesB: F64FlatArray,
        queryC: F64FlatArray,
        aggReducer: AggregationWithReducerWithArgWithCompare = MAX,
        lagBound: Int? = null,
        queryCIndex: Int = -1, // can be < 0 to make queryCIndex follow queryBIndex
    ): Double {
        val m = queryC.length
        val nB = timeSeriesB.length

        if (lagBound != null) {

            // the bound would make that no m values are left in TB
            if (queryCIndex > nB - m) {
                return aggReducer.initializer
            }

            val startBoundIndex = if (queryCIndex < 0) 0 else (queryCIndex - lagBound).coerceAtLeast(0)
            val endBoundIndex = if (queryCIndex < 0) nB - m else (queryCIndex + lagBound).coerceAtMost(nB - m)

//        val results = F64FlatArray.full()
            var best = aggReducer.initializer
            timeSeriesB
                .sliceFlat(from = startBoundIndex, to = endBoundIndex + m)
                .forEachWindowWithIndex(queryC.length) { queryBIndex, queryB ->
                    val result = pearson3dLine(
                        timeSeriesA = timeSeriesA,
                        queryB = queryB,
                        queryC = queryC,
                        aggReducer = aggReducer,
                        lagBound = lagBound,
                        queryBIndex = queryBIndex + startBoundIndex,
                        queryCIndex = if (queryCIndex < 0) queryBIndex else queryCIndex,
                    )
                    if (aggReducer.firstIsBetterThanSecond(first = result, second = best)) {
                        best = aggReducer(result, best)
                    }
                }

            return best
        } else {
            var best = aggReducer.initializer
            timeSeriesB.forEachWindow(queryC.length) { queryB ->
                val result = pearson3dLine(
                    timeSeriesA = timeSeriesA,
                    queryB = queryB,
                    queryC = queryC,
                    aggReducer = aggReducer
                )

                if (aggReducer.firstIsBetterThanSecond(first = result, second = best)) {
                    best = aggReducer(result, best)
                }
            }

            return best
        }
    }

    // 1 time series, 2 queries, lagbound
    private fun pearson3dLine(
        timeSeriesA: F64FlatArray,
        queryB: F64FlatArray,
        queryC: F64FlatArray,
        aggReducer: AggregationWithReducerWithArgWithCompare = MAX,
        lagBound: Int? = null,
        queryBIndex: Int = -1,
        queryCIndex: Int = -1,
    ): Double {
        val m = queryB.length
        val nA = timeSeriesA.length

        if (lagBound != null) {

            require(abs(queryBIndex - queryCIndex) <= lagBound) { "Illegal input: query1Index and query2Index are too far apart for the lagbound." }

            // the bound would make that no m values are left in TA
            if (queryCIndex > nA - m || queryBIndex > nA - m) {
                return aggReducer.initializer
            }

            val startBoundIndex = (maxOf(queryBIndex, queryCIndex) - lagBound).coerceAtLeast(0)
            val endBoundIndex = (minOf(queryBIndex, queryCIndex) + lagBound).coerceAtMost(nA - m)

            var best = aggReducer.initializer

            timeSeriesA
                .sliceFlat(from = startBoundIndex, to = endBoundIndex + m)
                .forEachWindow(queryB.length) { queryA ->

                    val a = Pearson.correlationFunction(
                        queryA.copy().also {
                            it += queryC
                            it /= 2.0
                        },
                        queryB,
                    )

                    val b = Pearson.correlationFunction(
                        queryA.copy().also {
                            it += queryB
                            it /= 2.0
                        },
                        queryC,
                    )

                    val c = Pearson.correlationFunction(
                        queryB.copy().also {
                            it += queryC
                            it /= 2.0
                        },
                        queryA,
                    )


                    val result = aggReducer(a, b, c).coerceNotNaN(aggReducer.initializer)
                    if (aggReducer.firstIsBetterThanSecond(first = result, second = best)) {
                        best = aggReducer(result, best)
                    }
                }

            return best
        } else {

            var best = aggReducer.initializer

            timeSeriesA.forEachWindow(queryB.length) { queryA ->

                val a = Pearson.correlationFunction(
                    queryA.copy().also {
                        it += queryC
                        it /= 2.0
                    },
                    queryB,
                )

                val b = Pearson.correlationFunction(
                    queryA.copy().also {
                        it += queryB
                        it /= 2.0
                    },
                    queryC,
                )

                val c = Pearson.correlationFunction(
                    queryB.copy().also {
                        it += queryC
                        it /= 2.0
                    },
                    queryA,
                )

                val result = aggReducer(a, b, c).coerceNotNaN(aggReducer.initializer)
                if (aggReducer.firstIsBetterThanSecond(first = result, second = best)) {
                    best = aggReducer(result, best)
                }
            }

            return best
        }
    }


}