@file:Suppress("LocalVariableName", "NAME_SHADOWING")

package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson


import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.AggregationWithReducerWithArgWithCompare
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampArrayCache.Slot
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor.asF64Array
import org.jtransforms.fft.DoubleFFT_1D
import java.io.Serializable
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * Implementation of STAMP-Pearson-3TS. Calculates the best lagged ternary correlation for three time series.
 *
 * @param arrayCache (optional) the array cache that is used for this [StampPearson3ts] instance.
 * @param fourierCache (optional) the [DoubleFFT_1D] instance cache used for this [StampPearson3ts] instance.
 */
class StampPearson3ts @JvmOverloads constructor(
    arrayCache: StampArrayCache = StampArrayCache(
        maxSize = 100,
        slotManager = StampArrayCache.SlotsWithReuse,
    ),
    fourierCache: MutableMap<Int, DoubleFFT_1D> = mutableMapOf(),
) : StampPearson(arrayCache, fourierCache), Serializable {

    /**
     * @param maxArraySize The size of the largest series put in the algorithms. The array cache automatically enlarges
     *  if not large enough.
     * @param slotManager The method of reusing slots in the array cache. Default is [StampArrayCache.SlotsWithReuse].
     *  Use [StampArrayCache.NoSlots] if [StampPearson3ts] instance is accessed across multiple threads.
     */
    @JvmOverloads
    constructor(
        maxArraySize: Int,
        slotManager: StampArrayCache.SlotManager = StampArrayCache.SlotsWithReuse,
    ) : this(
        arrayCache = StampArrayCache(
            maxSize = maxArraySize,
            slotManager = slotManager,
        ),
        fourierCache = mutableMapOf(),
    )

    enum class TimeSeries {
        A, B, C
    }

    enum class AggregationMethod(val textual: String, val aggregatedTimeSeries: List<TimeSeries>) :
        Serializable {
        AB("agg(a, b) and c", listOf(TimeSeries.A, TimeSeries.B)),
        BC("a and agg(b, c)", listOf(TimeSeries.B, TimeSeries.C)),
        AC("agg(a, c) and b", listOf(TimeSeries.A, TimeSeries.C)),
        NONE("Aggregation didn't happen", emptyList()),
    }

    /**
     * Find the best correlation between 3 time series for each possible window position in each.
     * ([DoubleArray] version)
     *
     * @param timeSeriesA First time series.
     * @param timeSeriesB Second time series.
     * @param timeSeriesC Third time series (will be split first).
     * @param windowSize, also called `m`, the window size used to traverse all three time series.
     * @param reducer Defines whether to find the highest or lowest correlation.
     *      In practice this is either [MIN] or [MAX]. Defaults to [MAX].
     *
     * @param timeSeriesASlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesA]. NOTE: As the name suggests, if [timeSeriesA] contains `NaN`s,
     *      use [timeSeriesAWithoutNaN] to calculate it.
     * @param timeSeriesASlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesA]. NOTE: As the name suggests, if [timeSeriesA] contains `NaN`s,
     *      use [timeSeriesAWithoutNaN] to calculate it.
     * @param timeSeriesBSlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesB]. It doesn't matter whether [timeSeriesB] has any `NaN`s.
     * @param timeSeriesBSlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesB]. It doesn't matter whether [timeSeriesB] has any `NaN`s.
     * @param timeSeriesCSlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesC]. It doesn't matter whether [timeSeriesC] has any `NaN`s.
     * @param timeSeriesCSlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesC]. It doesn't matter whether [timeSeriesC] has any `NaN`s.
     *
     * @param lagBound Optional enveloping mechanism. If provided, all window placements are within [lagBound] distance of each other.
     * @param resultCubeCallback Optional matrix that can be provided to be filled with all results.
     *      It's faster when not provided. Matrix must be of dimensions
     *      ([timeSeriesC]`.size` - [windowSize] + 1,
     *          [timeSeriesB]`.size` - [windowSize] + 1,
     *          [timeSeriesA]`.size` - [windowSize] + 1).
     * @param resultIndexCallback Optional callback function that can be provided which will be called before
     *      [massPearson3dCube] returns with the indices and aggregation method of the best found correlation.
     * @return The value of the best (according to [reducer]) correlation found between the three given time series.
     */
    @JvmOverloads
    fun stampPearson3ts(
        timeSeriesA: DoubleArray,
        timeSeriesB: DoubleArray,
        timeSeriesC: DoubleArray, // will be split first
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        timeSeriesASlidingMeans: DoubleArray? = null,
        timeSeriesASlidingStds: DoubleArray? = null,
        timeSeriesBSlidingMeans: DoubleArray? = null,
        timeSeriesBSlidingStds: DoubleArray? = null,
        timeSeriesCSlidingMeans: DoubleArray? = null,
        timeSeriesCSlidingStds: DoubleArray? = null,

        lagBound: Int? = null,

        // get optional result matrix
        resultCubeCallback: F64Array? = null,
        resultIndexCallback: ((aggregation: AggregationMethod, aIndex: Int, bIndex: Int, cIndex: Int) -> Unit)? = null,
    ): Double = stampPearson3ts(
        timeSeriesA = timeSeriesA.asF64Array(),
        timeSeriesB = timeSeriesB.asF64Array(),
        timeSeriesC = timeSeriesC.asF64Array(),
        windowSize = windowSize,
        reducer = reducer,
        timeSeriesASlidingMeans = timeSeriesASlidingMeans?.asF64Array(),
        timeSeriesASlidingStds = timeSeriesASlidingStds?.asF64Array(),
        timeSeriesBSlidingMeans = timeSeriesBSlidingMeans?.asF64Array(),
        timeSeriesBSlidingStds = timeSeriesBSlidingStds?.asF64Array(),
        timeSeriesCSlidingMeans = timeSeriesCSlidingMeans?.asF64Array(),
        timeSeriesCSlidingStds = timeSeriesCSlidingStds?.asF64Array(),
        lagBound = lagBound,
        resultCubeCallback = resultCubeCallback,
        resultIndexCallback = resultIndexCallback
    )

    /**
     * Find the best correlation between 3 time series for each possible window position in each.
     *
     * @param timeSeriesA First time series.
     * @param timeSeriesB Second time series.
     * @param timeSeriesC Third time series (will be split first).
     * @param windowSize, also called `m`, the window size used to traverse all three time series.
     * @param reducer Defines whether to find the highest or lowest correlation.
     *      In practice this is either [MIN] or [MAX]. Defaults to [MAX].
     *
     * @param timeSeriesASlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesA]. NOTE: As the name suggests, if [timeSeriesA] contains `NaN`s,
     *      use [timeSeriesAWithoutNaN] to calculate it.
     * @param timeSeriesASlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesA]. NOTE: As the name suggests, if [timeSeriesA] contains `NaN`s,
     *      use [timeSeriesAWithoutNaN] to calculate it.
     * @param timeSeriesBSlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesB]. It doesn't matter whether [timeSeriesB] has any `NaN`s.
     * @param timeSeriesBSlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesB]. It doesn't matter whether [timeSeriesB] has any `NaN`s.
     * @param timeSeriesCSlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesC]. It doesn't matter whether [timeSeriesC] has any `NaN`s.
     * @param timeSeriesCSlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesC]. It doesn't matter whether [timeSeriesC] has any `NaN`s.
     *
     * @param lagBound Optional enveloping mechanism. If provided, all window placements are within [lagBound] distance of each other.
     * @param resultCubeCallback Optional matrix that can be provided to be filled with all results.
     *      It's faster when not provided. Matrix must be of dimensions
     *      ([timeSeriesC]`.size` - [windowSize] + 1,
     *          [timeSeriesB]`.size` - [windowSize] + 1,
     *          [timeSeriesA]`.size` - [windowSize] + 1).
     * @param resultIndexCallback Optional callback function that can be provided which will be called before
     *      [massPearson3dCube] returns with the indices and aggregation method of the best found correlation.
     * @return The value of the best (according to [reducer]) correlation found between the three given time series.
     */
    @JvmOverloads
    fun stampPearson3ts(
        timeSeriesA: F64FlatArray,
        timeSeriesB: F64FlatArray,
        timeSeriesC: F64FlatArray, // will be split first
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        timeSeriesASlidingMeans: F64FlatArray? = null,
        timeSeriesASlidingStds: F64FlatArray? = null,
        timeSeriesBSlidingMeans: F64FlatArray? = null,
        timeSeriesBSlidingStds: F64FlatArray? = null,
        timeSeriesCSlidingMeans: F64FlatArray? = null,
        timeSeriesCSlidingStds: F64FlatArray? = null,

        lagBound: Int? = null,

        // get optional result matrix
        resultCubeCallback: F64Array? = null,
        resultIndexCallback: ((aggregation: AggregationMethod, aIndex: Int, bIndex: Int, cIndex: Int) -> Unit)? = null,
    ): Double {
        val m = windowSize

        val TA: F64FlatArray = timeSeriesA
        val TB: F64FlatArray = timeSeriesB
        val TC: F64FlatArray = timeSeriesC

        val nA: Int = TA.length
        val nB: Int = TB.length
        val nC: Int = TC.length


        if (resultCubeCallback != null) require(
            resultCubeCallback.shape.contentEquals(
                intArrayOf(
                    nC - m + 1,
                    nB - m + 1,
                    nA - m + 1,
                )
            )
        ) {
            "the result should be of size [nC - m + 1, nB - m + 1, nA - m + 1], where [n] is the size of [timeSeries] and [m] is the size of [query]."
        }

        // can be reused in lower levels

        var timeSeriesASlidingMeans = timeSeriesASlidingMeans
        var timeSeriesASlidingStds = timeSeriesASlidingStds

        if (timeSeriesASlidingStds == null || timeSeriesASlidingMeans == null) {
            timeSeriesASlidingMeans = arrayCache[nA - m + 1, Slot.F]
            timeSeriesASlidingStds = arrayCache[nA - m + 1, Slot.G]
            computeSlidingMeanStd(
                T = TA,
                windowSize = m,
                destMeanArray = timeSeriesASlidingMeans,
                destStdArray = timeSeriesASlidingStds,
            )
        }

        var timeSeriesBSlidingMeans = timeSeriesBSlidingMeans
        var timeSeriesBSlidingStds = timeSeriesBSlidingStds

        if (timeSeriesBSlidingMeans == null || timeSeriesBSlidingStds == null) {
            timeSeriesBSlidingMeans = arrayCache[nB - m + 1, Slot.J]
            timeSeriesBSlidingStds = arrayCache[nB - m + 1, Slot.K]
            computeSlidingMeanStd(
                T = TB,
                windowSize = m,
                destMeanArray = timeSeriesBSlidingMeans,
                destStdArray = timeSeriesBSlidingStds,
            )
        }

        var timeSeriesCSlidingMeans = timeSeriesCSlidingMeans
        var timeSeriesCSlidingStds = timeSeriesCSlidingStds

        if (timeSeriesCSlidingMeans == null || timeSeriesCSlidingStds == null) {
            timeSeriesCSlidingMeans = arrayCache[nC - m + 1, Slot.S]
            timeSeriesCSlidingStds = arrayCache[nC - m + 1, Slot.T]
            computeSlidingMeanStd(
                T = TC,
                windowSize = m,
                destMeanArray = timeSeriesCSlidingMeans,
                destStdArray = timeSeriesCSlidingStds,
            )
        }

        val considerCacheQueryCTimeSeriesASlidingDotProductsEmpty = delegateOf(true)
        val cacheQueryCTimeSeriesASlidingDotProducts = arrayCache[nA - m + 1, Slot.C]

        val considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty = delegateOf(true)
        val cacheQueryCTimeSeriesBSlidingDotProducts = arrayCache[nB - m + 1, Slot.D]

        val considerCacheQueryBTimeSeriesASlidingDotProductsEmpty = delegateOf(true)
        val cacheQueryBTimeSeriesASlidingDotProducts = arrayCache[nA - m + 1, Slot.E]

        var previousQueryC: F64FlatArray? = null

        var bestAggregation = AggregationMethod.NONE
        var bestAIndex = -1
        var bestBIndex = -1
        var bestCIndex = -1

        var bestCorrelation = reducer.initializer

        TC.forEachWindowWithIndex(windowSize = m) { queryCIndex, queryC ->

            // skip for continuity
            if (queryC.any { it.isNaN() }) {
                considerCacheQueryCTimeSeriesASlidingDotProductsEmpty.value = true
                considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty.value = true

                return@forEachWindowWithIndex
            }

            var aggregation = AggregationMethod.NONE
            var aIndex: Int = -1
            var bIndex: Int = -1

            val planeResult = stampPearson3tsSub(
                timeSeriesA = TA,
                timeSeriesB = TB,
                queryC = queryC,
                reducer = reducer,
                lagBound = lagBound,
                queryCIndex = queryCIndex,

                timeSeriesASlidingMeans = timeSeriesASlidingMeans,
                timeSeriesASlidingStds = timeSeriesASlidingStds,
                timeSeriesBSlidingMeans = timeSeriesBSlidingMeans,
                timeSeriesBSlidingStds = timeSeriesBSlidingStds,
                queryCMean = timeSeriesCSlidingMeans[queryCIndex],
                queryCStd = timeSeriesCSlidingStds[queryCIndex],
                considerCacheQueryCTimeSeriesASlidingDotProductsEmpty = considerCacheQueryCTimeSeriesASlidingDotProductsEmpty,
                cacheQueryCTimeSeriesASlidingDotProducts = cacheQueryCTimeSeriesASlidingDotProducts,
                considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty = considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty,
                cacheQueryCTimeSeriesBSlidingDotProducts = cacheQueryCTimeSeriesBSlidingDotProducts,
                considerCacheQueryBTimeSeriesASlidingDotProductsEmpty = considerCacheQueryBTimeSeriesASlidingDotProductsEmpty,
                cacheQueryBTimeSeriesASlidingDotProducts = cacheQueryBTimeSeriesASlidingDotProducts,
                previousQueryC = previousQueryC,

                resultPlaneCallback = resultCubeCallback?.view(index = queryCIndex),
                resultIndexCallback = when (resultIndexCallback) {
                    null -> null
                    else -> {
                        { aggI, aI, bI ->
                            aggregation = aggI
                            aIndex = aI
                            bIndex = bI
                        }
                    }
                },
            )

            if (reducer.firstIsBetterThanSecond(first = planeResult, second = bestCorrelation)) {
                bestCorrelation = reducer(bestCorrelation, planeResult)

                @Suppress("SENSELESS_COMPARISON")
                if (resultIndexCallback != null) { // don't listen to this hint, it's wrong
                    bestAggregation = aggregation
                    bestAIndex = aIndex
                    bestBIndex = bIndex
                    bestCIndex = queryCIndex
                }
            }
            previousQueryC = queryC
        }

        resultIndexCallback?.invoke(bestAggregation, bestAIndex, bestBIndex, bestCIndex)

        return bestCorrelation
    }

    internal fun stampPearson3tsSub(
        timeSeriesA: F64FlatArray,
        timeSeriesB: F64FlatArray,  // will be split up in queries
        queryC: F64FlatArray,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        // optional data for the lag bound
        lagBound: Int? = null,
        queryCIndex: Int? = null,

        // optional vertical optimizations
        timeSeriesASlidingMeans: F64FlatArray? = null,
        timeSeriesASlidingStds: F64FlatArray? = null,
        timeSeriesBSlidingMeans: F64FlatArray? = null,
        timeSeriesBSlidingStds: F64FlatArray? = null,
        queryCMean: Double? = null,
        queryCStd: Double? = null,

        // optional horizontal optimizations
        considerCacheQueryCTimeSeriesASlidingDotProductsEmpty: Delegate<Boolean>? = null,
        cacheQueryCTimeSeriesASlidingDotProducts: F64FlatArray? = null,
        considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty: Delegate<Boolean>? = null,
        cacheQueryCTimeSeriesBSlidingDotProducts: F64FlatArray? = null,
        considerCacheQueryBTimeSeriesASlidingDotProductsEmpty: Delegate<Boolean>? = null,
        cacheQueryBTimeSeriesASlidingDotProducts: F64FlatArray? = null,
        previousQueryC: F64FlatArray? = null,

        // get optional result plane
        resultPlaneCallback: F64Array? = null,
        resultIndexCallback: ((aggregation: AggregationMethod, aIndex: Int, bIndex: Int) -> Unit)? = null,
    ): Double {
        val TA = timeSeriesA
        var TB = timeSeriesB

        val m = queryC.length
        val nA = TA.length
        var nB = TB.length

        if (resultPlaneCallback != null) require(
            resultPlaneCallback.shape.contentEquals(
                intArrayOf(
                    nB - m + 1,
                    nA - m + 1,
                )
            )
        ) {
            "the result should be of size [nB - m + 1, nA - m + 1], where [n] is the size of [timeSeries] and [m] is the size of [query].\n" +
                    "currently it's ${resultPlaneCallback.shape.toList()}, nA=$nA, nB=$nB, m=$m"
        }

        // make input relevant to TB variable and apply lag bound to all precalculated values
        var timeSeriesBSlidingMeans = timeSeriesBSlidingMeans
        var timeSeriesBSlidingStds = timeSeriesBSlidingStds
        var cacheQueryCTimeSeriesBSlidingDotProducts = cacheQueryCTimeSeriesBSlidingDotProducts
        var bIndexModifier = 0

        if (lagBound != null) {

            // the bound would make that no m values are left in TB
            if (queryCIndex!! > nB - m) return reducer.initializer

            val startBoundIndex = (queryCIndex - lagBound).coerceAtLeast(0)
            val endBoundIndex = (queryCIndex + lagBound).coerceAtMost(nB - m)

            TB = TB.sliceFlat(from = startBoundIndex, to = endBoundIndex + m)

            nB = TB.length
            bIndexModifier = startBoundIndex

            timeSeriesBSlidingMeans = timeSeriesBSlidingMeans
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            timeSeriesBSlidingStds = timeSeriesBSlidingStds
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            cacheQueryCTimeSeriesBSlidingDotProducts = cacheQueryCTimeSeriesBSlidingDotProducts
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)
        }

        // reused later in agg(T, Q1) and Q2
        //      where Q1 is from [timeSeriesToBeSplit] and Q2 is [query]
        val queryCTimeSeriesASlidingDotProducts = when {
            cacheQueryCTimeSeriesASlidingDotProducts != null && considerCacheQueryCTimeSeriesASlidingDotProductsEmpty != null ->
                cacheQueryCTimeSeriesASlidingDotProducts.also {
                    updateQTCache(
                        QTCache = it,
                        considerQTCacheEmpty = considerCacheQueryCTimeSeriesASlidingDotProductsEmpty,
                        Q = queryC,
                        T = TA,
                        n = nA,
                        m = m,
                        previousQuery = previousQueryC,
                        returnArraySlot = Slot.B1,
                    )
                }
            else -> slidingDotProducts(Q = queryC, T = TA, returnArraySlot = Slot.B1)
        }

        // precalculate all QB dot QC
        val queryCTimeSeriesBSlidingDotProducts = when {
            cacheQueryCTimeSeriesBSlidingDotProducts != null && considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty != null ->
                cacheQueryCTimeSeriesBSlidingDotProducts.also {
                    updateQTCache(
                        QTCache = it,
                        considerQTCacheEmpty = considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty,
                        Q = queryC,
                        T = TB,
                        n = nB,
                        m = m,
                        previousQuery = previousQueryC,
                        returnArraySlot = Slot.B2,
                    )
                }
            else -> slidingDotProducts(Q = queryC, T = TB, returnArraySlot = Slot.B2)
        }

        // reused later in agg(T, Q1) and Q2
        //      where Q1 is from [timeSeriesToBeSplit] and Q2 is [query]
        val queryCMean = queryCMean ?: queryC.mean()
        val queryCStd = queryCStd ?: arrayCache[queryC.length, Slot.M] { queryC[it] }
            .let {
                it -= queryCMean
                sqrt(it.dot() / it.length)
            }

        // reused later in T and agg(Q1, Q2)
        //      where Q1 is from [timeSeriesToBeSplit] and Q2 is [query]
        var timeSeriesASlidingMeans = timeSeriesASlidingMeans
        var timeSeriesASlidingStds = timeSeriesASlidingStds

        if (timeSeriesASlidingStds == null || timeSeriesASlidingMeans == null) {
            timeSeriesASlidingMeans = arrayCache[nA - m + 1, Slot.F]
            timeSeriesASlidingStds = arrayCache[nA - m + 1, Slot.G]
            computeSlidingMeanStd(
                T = TA, // no difference with TA in practice
                windowSize = m,
                destMeanArray = timeSeriesASlidingMeans,
                destStdArray = timeSeriesASlidingStds,
            )
        }

        val timeSeriesASlidingMeansAggQC = arrayCache[nA - m + 1, Slot.H]
        val timeSeriesASlidingStdsAggQC = arrayCache[nA - m + 1, Slot.I]
        computeSlidingMeanStdWithAgg(
            T = TA,
            windowSize = m,
            destMeanArray = timeSeriesASlidingMeansAggQC,
            destStdArray = timeSeriesASlidingStdsAggQC,
            precalculatedTMeans = timeSeriesASlidingMeans,
            aggregations = queryC,
        )

        if (timeSeriesBSlidingMeans == null || timeSeriesBSlidingStds == null) {
            timeSeriesBSlidingMeans = arrayCache[nB - m + 1, Slot.J]
            timeSeriesBSlidingStds = arrayCache[nB - m + 1, Slot.K]
            computeSlidingMeanStd(
                T = TB,
                windowSize = m,
                destMeanArray = timeSeriesBSlidingMeans,
                destStdArray = timeSeriesBSlidingStds,
            )
        }

        val considerCacheQueryBTimeSeriesASlidingDotProductsEmpty =
            considerCacheQueryBTimeSeriesASlidingDotProductsEmpty ?: delegateOf(true)
        val cacheQueryBTimeSeriesASlidingDotProducts =
            if (cacheQueryBTimeSeriesASlidingDotProducts != null) {
                considerCacheQueryBTimeSeriesASlidingDotProductsEmpty.value = true
                cacheQueryBTimeSeriesASlidingDotProducts
            } else {
                arrayCache[nA - m + 1, Slot.L]
            }


        var bestAggregation: AggregationMethod = AggregationMethod.NONE
        var bestAIndex: Int = -1
        var bestBIndex: Int = -1

        var bestCorrelation: Double = reducer.initializer
        var previousQueryB: F64FlatArray? = null
        TB.forEachWindowWithIndex(m) { queryBIndex, queryB ->

            // skip for continuity
            if (queryB.any { it.isNaN() }) {
                considerCacheQueryBTimeSeriesASlidingDotProductsEmpty.value = true

                return@forEachWindowWithIndex
            }

            var aggregation: AggregationMethod = AggregationMethod.NONE
            var aIndex: Int = -1

            val lineResult = massPearson3ts(
                timeSeriesA = TA,
                queryB = queryB,
                queryC = queryC,
                reducer = reducer,

                lagBound = lagBound,
                queryBIndex = queryBIndex + bIndexModifier,
                queryCIndex = queryCIndex,

                queryCTimeSeriesASlidingDotProducts = queryCTimeSeriesASlidingDotProducts,
                queryBMean = timeSeriesBSlidingMeans[queryBIndex],
                queryBStd = timeSeriesBSlidingStds[queryBIndex],
                queryCMean = queryCMean,
                queryCStd = queryCStd,
                timeSeriesASlidingMeans = timeSeriesASlidingMeans,
                timeSeriesASlidingStds = timeSeriesASlidingStds,
                timeSeriesASlidingMeansAggQC = timeSeriesASlidingMeansAggQC,
                timeSeriesASlidingStdsAggQC = timeSeriesASlidingStdsAggQC,
                queryBDotQueryC = queryCTimeSeriesBSlidingDotProducts[queryBIndex],
                considerCacheQueryBTimeSeriesASlidingDotProductsEmpty = considerCacheQueryBTimeSeriesASlidingDotProductsEmpty,
                cacheQueryBTimeSeriesASlidingDotProducts = cacheQueryBTimeSeriesASlidingDotProducts,
                previousQueryB = previousQueryB,

                resultLineCallback = resultPlaneCallback?.viewAsFlat(queryBIndex + bIndexModifier),
                resultIndexCallback = when (resultIndexCallback) {
                    null -> null
                    else -> {
                        { aggI, aI ->
                            aggregation = aggI
                            aIndex = aI
                        }
                    }
                },
            )

            if (reducer.firstIsBetterThanSecond(first = lineResult, second = bestCorrelation)) {
                bestCorrelation = reducer(bestCorrelation, lineResult)

                @Suppress("SENSELESS_COMPARISON")
                if (resultIndexCallback != null) { // don't listen to this hint, it's wrong
                    bestAggregation = aggregation
                    bestAIndex = aIndex
                    bestBIndex = queryBIndex
                }
            }
            previousQueryB = queryB
        }

        resultIndexCallback?.invoke(bestAggregation, bestAIndex, bestBIndex + bIndexModifier)

        return bestCorrelation
    }

    @Suppress("SENSELESS_COMPARISON")
    internal fun massPearson3ts(
        timeSeriesA: F64FlatArray,
        queryB: F64FlatArray,
        queryC: F64FlatArray,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        // optional data for the lag bound
        lagBound: Int? = null,
        queryBIndex: Int? = null,
        queryCIndex: Int? = null,

        // optional vertical optimizations
        queryCTimeSeriesASlidingDotProducts: F64FlatArray? = null,
        queryBMean: Double? = null,
        queryBStd: Double? = null,
        queryCMean: Double? = null,
        queryCStd: Double? = null,
        timeSeriesASlidingMeans: F64FlatArray? = null,
        timeSeriesASlidingStds: F64FlatArray? = null,
        timeSeriesASlidingMeansAggQC: F64FlatArray? = null,
        timeSeriesASlidingStdsAggQC: F64FlatArray? = null,
        queryBDotQueryC: Double? = null,

        // optional horizontal optimizations
        considerCacheQueryBTimeSeriesASlidingDotProductsEmpty: Delegate<Boolean>? = null,
        cacheQueryBTimeSeriesASlidingDotProducts: F64FlatArray? = null,
        previousQueryB: F64FlatArray? = null,

        // get optional result line
        resultLineCallback: F64FlatArray? = null,

        // get optional result index
        resultIndexCallback: ((aggregation: AggregationMethod, aIndex: Int) -> Unit)? = null,
    ): Double {
        val QB = queryB
        val QC = queryC

        require(QB.length == QC.length)

        val m = QB.length

        var TA = timeSeriesA
        var nA = TA.length

        require(TA.all { it.isNotNaN() }) { "timeSeriesA is not allowed to have continuity removed and thus may not contain NaNs." }

        if (resultLineCallback != null) require(resultLineCallback.length == nA - m + 1) {
            "the result should be of size [nA - m + 1], where [n] is the size of [timeSeries] and [m] is the size of [query]."
        }

        // make input relevant to TA variable and apply lag bound to all precalculated values
        var timeSeriesASlidingMeans = timeSeriesASlidingMeans
        var timeSeriesASlidingStds = timeSeriesASlidingStds
        var timeSeriesASlidingMeansAggQC = timeSeriesASlidingMeansAggQC
        var timeSeriesASlidingStdsAggQC = timeSeriesASlidingStdsAggQC
        var cacheQueryBTimeSeriesASlidingDotProducts = cacheQueryBTimeSeriesASlidingDotProducts
        var queryCTimeSeriesASlidingDotProducts = queryCTimeSeriesASlidingDotProducts
        var aIndexModifier = 0

        if (lagBound != null) {
            require(abs(queryBIndex!! - queryCIndex!!) <= lagBound) { "Illegal input: queryBIndex and queryCIndex are too far apart for the lagbound." }

            // the bound would make that no m values are left in TA
            if (queryCIndex > nA - m || queryBIndex > nA - m) return reducer.initializer

            val startBoundIndex = (maxOf(queryBIndex, queryCIndex) - lagBound).coerceAtLeast(0)
            val endBoundIndex = (minOf(queryBIndex, queryCIndex) + lagBound).coerceAtMost(nA - m)

            TA = TA.sliceFlat(from = startBoundIndex, to = endBoundIndex + m)
            nA = TA.length
            aIndexModifier = startBoundIndex

            timeSeriesASlidingMeans = timeSeriesASlidingMeans
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            timeSeriesASlidingStds = timeSeriesASlidingStds
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            timeSeriesASlidingMeansAggQC = timeSeriesASlidingMeansAggQC
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            timeSeriesASlidingStdsAggQC = timeSeriesASlidingStdsAggQC
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            queryCTimeSeriesASlidingDotProducts = queryCTimeSeriesASlidingDotProducts
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)

            cacheQueryBTimeSeriesASlidingDotProducts = cacheQueryBTimeSeriesASlidingDotProducts
                ?.sliceFlat(from = startBoundIndex, to = endBoundIndex + 1)
        }


        // used in  agg(TA, QB) and QC   and   agg(TA, QC) and QB
        val QBdotQC = queryBDotQueryC ?: (QB dot QC)

        val resultsMatrix: F64Array? = if (resultLineCallback == null) null else F64Array(3, nA - m + 1)
        val resultsIndexArray: IntArray? = if (resultIndexCallback == null) null else IntArray(3) { -1 }

        // agg(TA, QB) and QC
        val a = massPearson(
            timeSeries = TA,
            aggregationQueries = QB,
            query = arrayCache[QC.length, Slot.U].also {
                QC.copyTo(it)
            },
            reducer = reducer,

            queryDotAggregationProducts = QBdotQC,
            queryTimeSeriesSlidingDotProducts = queryCTimeSeriesASlidingDotProducts,
            queryMean = queryCMean,
            queryStd = queryCStd,
            precalculatedTimeSeriesSlidingMeans = timeSeriesASlidingMeans,

            resultLineCallback = resultsMatrix?.let { _ ->
                fun(resultIndex: F64FlatArray) {
                    resultsMatrix.V[0] = resultIndex
                }
            },

            resultIndexCallback = resultsIndexArray?.let { _ ->
                fun(resultIndex: Int) {
                    resultsIndexArray[0] = resultIndex
                }
            },
        ).coerceNotNaN(reducer.initializer)

        // TA and agg(QB, QC) NOTE, is done without aggregation, but can use the same function
        val b = massPearson(
            timeSeries = TA,
            query = arrayCache[QB.length, Slot.N].also {
                QB.copyTo(it)
                it += QC
                it /= 2.0
            },
            reducer = reducer,

            timeSeriesSlidingMeans = timeSeriesASlidingMeans,
            timeSeriesSlidingStds = timeSeriesASlidingStds,

            resultLineCallback = resultsMatrix?.let { _ ->
                fun(resultIndex: F64FlatArray) {
                    resultsMatrix.V[1] = resultIndex
                }
            },
            resultIndexCallback = resultsIndexArray?.let { _ ->
                fun(resultIndex: Int) {
                    resultsIndexArray[1] = resultIndex
                }
            },
        ).coerceNotNaN(reducer.initializer)

        // agg(TA, QC) and QB
        val c = massPearson(
            timeSeries = TA,
            aggregationQueries = QC,
            query = QB,
            reducer = reducer,

            queryDotAggregationProducts = QBdotQC,
            timeSeriesSlidingMeans = timeSeriesASlidingMeansAggQC,
            timeSeriesSlidingStds = timeSeriesASlidingStdsAggQC,
            considerCacheQueryTimeSeriesSlidingDotProductsEmpty = considerCacheQueryBTimeSeriesASlidingDotProductsEmpty,
            cacheQueryTimeSeriesSlidingDotProducts = cacheQueryBTimeSeriesASlidingDotProducts,
            previousQuery = previousQueryB,
            queryMean = queryBMean,
            queryStd = queryBStd,

            resultLineCallback = resultsMatrix?.let { _ ->
                fun(resultIndex: F64FlatArray) {
                    resultsMatrix.V[2] = resultIndex
                }
            },
            resultIndexCallback = resultsIndexArray?.let { _ ->
                fun(resultIndex: Int) {
                    resultsIndexArray[2] = resultIndex
                }
            },
        ).coerceNotNaN(reducer.initializer)


        val bestCorrelation = reducer(a, b, c)

        if (resultIndexCallback != null) {
            val aggregation = when (bestCorrelation) {
                a -> 0
                b -> 1
                c -> 2
                else -> throw Exception("broken reducer? a: $a, b: $b, c: $c, reduced: $bestCorrelation")
            }
            val aIndex: Int = resultsIndexArray!![aggregation]

            resultIndexCallback(AggregationMethod.values()[aggregation], aIndex + aIndexModifier)
        }


        if (resultLineCallback != null) {
            resultsMatrix!!.transformInPlace { it.coerceNotNaN(reducer.initializer) }
            var i = 0
            resultsMatrix.viewAsFlat(0, 0).transformInPlace {
                reducer(
                    resultsMatrix.viewAsFlat(i++, 1)
                )
            }

            val line = resultsMatrix.viewAsFlat(0, 0)
            line.copyTo(resultLineCallback.sliceFlat(from = aIndexModifier, to = aIndexModifier + line.length))
        }

        return bestCorrelation
    }
}