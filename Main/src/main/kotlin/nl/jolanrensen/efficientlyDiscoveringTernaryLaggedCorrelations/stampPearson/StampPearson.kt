@file:Suppress("NOTHING_TO_INLINE", "NAME_SHADOWING")

package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampArrayCache.Slot
import org.jetbrains.bio.viktor.*
import org.jtransforms.fft.DoubleFFT_1D
import java.io.Serializable
import kotlin.math.sqrt

/**
 * Implementation of STAMP from Matrix Profile paper using Mueen's similarity search.
 * Calculates Pearson correlation instead of Euclidean distance, but can be easily converted using [StampPearson.pearsonCorrelationsToEuclidean].
 *
 * @param arrayCache (optional) the array cache that is used for this [StampPearson] instance.
 * @param fourierCache (optional) the [DoubleFFT_1D] instance cache used for this [StampPearson] instance.
 */
@Suppress("LocalVariableName")
open class StampPearson @JvmOverloads constructor(
    internal val arrayCache: StampArrayCache = StampArrayCache(
        maxSize = 100,
        slotManager = StampArrayCache.SlotsWithReuse,
    ),
    internal val fourierCache: MutableMap<Int, DoubleFFT_1D> = mutableMapOf(),
) : Serializable {

    /**
     * @param maxArraySize The size of the largest series put in the algorithms. The array cache automatically enlarges
     *  if not large enough.
     * @param slotManager The method of reusing slots in the array cache. Default is [StampArrayCache.SlotsWithReuse].
     *  Use [StampArrayCache.NoSlots] if [StampPearson] instance is accessed across multiple threads.
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


    /**
     * The Matrix Profile function using Pearson correlation ([DoubleArray] variant). It calculates the Matrix Profile for [timeSeries] in Pearson correlations.
     * It supports having straight-line skipping and automatically re-orders A and B for efficiency.
     * @see [StampPearson.pearsonCorrelationsToEuclidean] to get the Euclidean distances instead.
     *
     * @param timeSeries the time series. If skipping occurs in the series, also supply [timeSeriesWithoutNaN].
     * @param timeSeriesWithoutNaN (optional) the original time series that needs to be supplied if skipping occurs in [timeSeries].
     * @param windowSize the window size.
     * @param reducer (optional, default = MAX) either [MAX] or [MIN].
     *
     * @param timeSeriesSlidingMeans (optional) the sliding means of [timeSeries], so that it doesn't have to be
     *  calculated in the function anymore.
     * @param timeSeriesSlidingStds (optional) the sliding standard deviations of [timeSeries], so that it doesn't have
     *  to be calculated in the function anymore.
     *
     * @return the Pearson Matrix Profile and indices of the Matrix Profile for [timeSeries].
     */
    @JvmOverloads
    fun matrixProfilePearson(
        timeSeries: DoubleArray,
        timeSeriesWithoutNaN: DoubleArray? = null, // if any continuity is removed from [timeSeries], this needs to be supplied
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        // optional horizontal speed optimization
        timeSeriesSlidingMeans: F64FlatArray? = null,
        timeSeriesSlidingStds: F64FlatArray? = null,
    ): Pair<F64FlatArray, IntArray> =
        matrixProfilePearson(
            timeSeries = timeSeries.asF64Array(),
            timeSeriesWithoutNaN = timeSeriesWithoutNaN?.asF64Array(),
            windowSize = windowSize,
            reducer = reducer,
            timeSeriesSlidingMeans = timeSeriesSlidingMeans,
            timeSeriesSlidingStds = timeSeriesSlidingStds,
        )

    /**
     * The Matrix Profile function using Pearson correlation. It calculates the Matrix Profile for [timeSeries] in Pearson correlations.
     * It supports having straight-line skipping and automatically re-orders A and B for efficiency.
     * @see [StampPearson.pearsonCorrelationsToEuclidean] to get the Euclidean distances instead.
     *
     * @param timeSeries the time series. If skipping occurs in the series, also supply [timeSeriesWithoutNaN].
     * @param timeSeriesWithoutNaN (optional) the original time series that needs to be supplied if skipping occurs in [timeSeries].
     * @param windowSize the window size.
     * @param reducer (optional, default = MAX) either [MAX] or [MIN].
     *
     * @param timeSeriesSlidingMeans (optional) the sliding means of [timeSeries], so that it doesn't have to be
     *  calculated in the function anymore.
     * @param timeSeriesSlidingStds (optional) the sliding standard deviations of [timeSeries], so that it doesn't have
     *  to be calculated in the function anymore.
     *
     * @return the Pearson Matrix Profile and indices of the Matrix Profile for [timeSeries].
     */
    @JvmOverloads
    fun matrixProfilePearson(
        timeSeries: F64FlatArray,
        timeSeriesWithoutNaN: F64FlatArray? = null, // if any continuity is removed from [timeSeries], this needs to be supplied
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        // optional horizontal speed optimization
        timeSeriesSlidingMeans: F64FlatArray? = null,
        timeSeriesSlidingStds: F64FlatArray? = null,
    ): Pair<F64FlatArray, IntArray> {
        val n = timeSeries.length
        val m = windowSize

        if (timeSeries.any { it.isNaN() }) require(timeSeriesWithoutNaN != null) { "The supplied time series has continuity removed, please supply the original time series without NaNs as well." }

        var timeSeriesSlidingMeans = timeSeriesSlidingMeans
        var timeSeriesSlidingStds = timeSeriesSlidingStds

        if (timeSeriesSlidingMeans == null || timeSeriesSlidingStds == null) {
            timeSeriesSlidingMeans = F64FlatArray(n - m + 1)
            timeSeriesSlidingStds = F64FlatArray(n - m + 1)
            computeSlidingMeanStd(
                T = timeSeries,
                windowSize = m,
                destMeanArray = timeSeriesSlidingMeans,
                destStdArray = timeSeriesSlidingStds,
            )
        }

        val planeResult = try {
            F64Array.full(
                n - m + 1,
                n - m + 1,
                init = reducer.initializer,
            )
        } catch (e: OutOfMemoryError) {
            throw Error("trying to create ${n - m + 1}^2 array", e)
        }
        stampPearson(
            timeSeriesA = timeSeries.copy(),
            timeSeriesB = timeSeries.copy(),
            windowSize = windowSize,

            timeSeriesAWithoutNaN = timeSeriesWithoutNaN?.copy(),
            timeSeriesBWithoutNaN = timeSeriesWithoutNaN?.copy(),
            reducer = reducer,

            timeSeriesASlidingMeans = timeSeriesSlidingMeans,
            timeSeriesASlidingStds = timeSeriesSlidingStds,
            timeSeriesBSlidingMeans = timeSeriesSlidingMeans,
            timeSeriesBSlidingStds = timeSeriesSlidingStds,

            resultPlaneCallback = planeResult,
        )

        // obviously the diagonal is not needed
        for (i in 0..n - m) {
            planeResult[i, i] = reducer.initializer
        }
        // fix infs
        planeResult.transformInPlace { it.coerceNotInf(reducer.initializer) }

        val lineResult = F64FlatArray(n - m + 1)
        val indexResult = IntArray(n - m + 1)

        for (i in 0..n - m) {
            val line = planeResult.viewAsFlat(i)
            val bestIndex = reducer.argFunction(line)
            indexResult[i] = bestIndex
            lineResult[i] = reducer(line)
        }

        return lineResult to indexResult
    }

    /**
     * The STAMP-Pearson function ([DoubleArray] variant). It finds the best Pearson correlation value between all window position in [timeSeriesA]
     * and [timeSeriesB]. It supports having straight-line skipping and automatically re-orders A and B for efficiency.
     *
     * @param timeSeriesA the first time series. If skipping occurs in the series, also supply [timeSeriesAWithoutNaN].
     * @param timeSeriesB the second time series. If skipping occurs in the series, also supply [timeSeriesBWithoutNaN].
     * @param windowSize the window size.
     * @param reducer (optional, default = MAX) either [MAX] or [MIN].
     * @param timeSeriesAWithoutNaN (optional) the original first time series that needs to be supplied if skipping occurs in [timeSeriesA].
     * @param timeSeriesBWithoutNaN (optional) the original second time series that needs to be supplied if skipping occurs in [timeSeriesB].
     *
     * @param timeSeriesASlidingMeans (optional) the sliding means of [timeSeriesA], so that it doesn't have to be
     *  calculated in the function anymore.
     * @param timeSeriesASlidingStds (optional) the sliding standard deviations of [timeSeriesA], so that it doesn't have
     *  to be calculated in the function anymore.
     * @param timeSeriesBSlidingMeans (optional) the sliding means of [timeSeriesB], so that it doesn't have to be
     *  calculated in the function anymore.
     * @param timeSeriesBSlidingStds (optional) the sliding standard deviations of [timeSeriesB], so that it doesn't have
     *  to be calculated in the function anymore.
     *
     * @param resultPlaneCallback (optional) this 2D array can be supplied to store all the results instead of just the best.
     *  Should be of size [nB - m + 1, nA - m + 1].
     * @param resultIndexCallback (optional) a callback to get the index of the start of the window positions in [timeSeriesA] and [timeSeriesB] for the best found correlation.
     *
     * @return the value of the best Pearson correlation.
     */
    @JvmOverloads
    fun stampPearson(
        timeSeriesA: DoubleArray,
        timeSeriesB: DoubleArray,
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,
        timeSeriesAWithoutNaN: DoubleArray? = null, // if any continuity is removed from [timeSeriesA], this needs to be supplied
        timeSeriesBWithoutNaN: DoubleArray? = null, // if any continuity is removed from [timeSeriesB], this needs to be supplied

        // optional horizontal optimizations
        timeSeriesASlidingMeans: DoubleArray? = null,
        timeSeriesASlidingStds: DoubleArray? = null,
        timeSeriesBSlidingMeans: DoubleArray? = null,
        timeSeriesBSlidingStds: DoubleArray? = null,

        // get optional result plane
        resultPlaneCallback: F64Array? = null,

        // get optional result index
        resultIndexCallback: ((aIndex: Int, bIndex: Int) -> Unit)? = null,
    ): Double = stampPearson(
        timeSeriesA = timeSeriesA.asF64Array(),
        timeSeriesB = timeSeriesB.asF64Array(),
        windowSize = windowSize,
        reducer = reducer,
        timeSeriesAWithoutNaN = timeSeriesAWithoutNaN?.asF64Array(),
        timeSeriesBWithoutNaN = timeSeriesBWithoutNaN?.asF64Array(),
        timeSeriesASlidingMeans = timeSeriesASlidingMeans?.asF64Array(),
        timeSeriesASlidingStds = timeSeriesASlidingStds?.asF64Array(),
        timeSeriesBSlidingMeans = timeSeriesBSlidingMeans?.asF64Array(),
        timeSeriesBSlidingStds = timeSeriesBSlidingStds?.asF64Array(),
        resultPlaneCallback = resultPlaneCallback,
        resultIndexCallback = resultIndexCallback,
    )

    /**
     * The STAMP-Pearson function. It finds the best Pearson correlation value between all window position in [timeSeriesA]
     * and [timeSeriesB]. It supports having straight-line skipping and automatically re-orders A and B for efficiency.
     *
     * @param timeSeriesA the first time series. If skipping occurs in the series, also supply [timeSeriesAWithoutNaN].
     * @param timeSeriesB the second time series. If skipping occurs in the series, also supply [timeSeriesBWithoutNaN].
     * @param windowSize the window size.
     * @param reducer (optional, default = MAX) either [MAX] or [MIN].
     * @param timeSeriesAWithoutNaN (optional) the original first time series that needs to be supplied if skipping occurs in [timeSeriesA].
     * @param timeSeriesBWithoutNaN (optional) the original second time series that needs to be supplied if skipping occurs in [timeSeriesB].
     *
     * @param timeSeriesASlidingMeans (optional) the sliding means of [timeSeriesA], so that it doesn't have to be
     *  calculated in the function anymore.
     * @param timeSeriesASlidingStds (optional) the sliding standard deviations of [timeSeriesA], so that it doesn't have
     *  to be calculated in the function anymore.
     * @param timeSeriesBSlidingMeans (optional) the sliding means of [timeSeriesB], so that it doesn't have to be
     *  calculated in the function anymore.
     * @param timeSeriesBSlidingStds (optional) the sliding standard deviations of [timeSeriesB], so that it doesn't have
     *  to be calculated in the function anymore.
     *
     * @param resultPlaneCallback (optional) this 2D array can be supplied to store all the results instead of just the best.
     *  Should be of size [nB - m + 1, nA - m + 1].
     * @param resultIndexCallback (optional) a callback to get the index of the start of the window positions in [timeSeriesA] and [timeSeriesB] for the best found correlation.
     *
     * @return the value of the best Pearson correlation.
     */
    @JvmOverloads
    fun stampPearson(
        timeSeriesA: F64FlatArray,
        timeSeriesB: F64FlatArray,
        windowSize: Int,
        reducer: AggregationWithReducerWithArgWithCompare = MAX,
        timeSeriesAWithoutNaN: F64FlatArray? = null, // if any continuity is removed from [timeSeriesA], this needs to be supplied
        timeSeriesBWithoutNaN: F64FlatArray? = null, // if any continuity is removed from [timeSeriesB], this needs to be supplied

        // optional horizontal optimizations
        timeSeriesASlidingMeans: F64FlatArray? = null,
        timeSeriesASlidingStds: F64FlatArray? = null,
        timeSeriesBSlidingMeans: F64FlatArray? = null,
        timeSeriesBSlidingStds: F64FlatArray? = null,

        // get optional result plane
        resultPlaneCallback: F64Array? = null,

        // get optional result index
        resultIndexCallback: ((aIndex: Int, bIndex: Int) -> Unit)? = null,
    ): Double {
        var TA = timeSeriesA
        var TB = timeSeriesB

        val TAHasNaNs = TA.any { it.isNaN() }
        val TBHasNaNs = TB.any { it.isNaN() }

        if (TAHasNaNs) require(timeSeriesAWithoutNaN != null) { "The supplied time series A has continuity removed, please supply the original time series without NaNs as well." }
        if (TBHasNaNs) require(timeSeriesBWithoutNaN != null) { "The supplied time series B has continuity removed, please supply the original time series without NaNs as well." }

        var TANoNaN = timeSeriesAWithoutNaN ?: TA
        var TBNoNaN = timeSeriesBWithoutNaN ?: TB

        // reorder A and B according to the amount of NaNs
        var resultIndexCallback = resultIndexCallback
        var aAndBAreFlipped = false
        if ((TAHasNaNs || TBHasNaNs) && TB.count { it.isNaN() } > TA.count { it.isNaN() }) {
            aAndBAreFlipped = true

            val tempT = TA
            TA = TB
            TB = tempT

            val tempTNoNaN = TANoNaN
            TANoNaN = TBNoNaN
            TBNoNaN = tempTNoNaN

            if (resultIndexCallback != null) {
                val oldResultIndexCallback = resultIndexCallback
                resultIndexCallback = { a, b ->
                    oldResultIndexCallback(b, a)
                }
            }
        }

        val nA = TA.length
        val nB = TB.length
        val m = windowSize

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

        var timeSeriesASlidingMeans = timeSeriesASlidingMeans
        var timeSeriesASlidingStds = timeSeriesASlidingStds

        if (timeSeriesASlidingStds == null || timeSeriesASlidingMeans == null) {
            timeSeriesASlidingMeans = F64FlatArray(nA - m + 1)
            timeSeriesASlidingStds = F64FlatArray(nA - m + 1)
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
            timeSeriesBSlidingMeans = F64FlatArray(nB - m + 1)
            timeSeriesBSlidingStds = F64FlatArray(nB - m + 1)
            computeSlidingMeanStd(
                T = TB,
                windowSize = m,
                destMeanArray = timeSeriesBSlidingMeans,
                destStdArray = timeSeriesBSlidingStds,
            )
        }

        val considerCacheQueryBTimeSeriesASlidingDotProductsEmpty = delegateOf(true)
        val cacheQueryBTimeSeriesASlidingDotProducts = F64FlatArray(nA - m + 1)

        fun getResultLineCallback(queryBIndex: Int): ((F64FlatArray) -> Unit)? =
            when (resultPlaneCallback) {
                null -> null
                else -> { it: F64FlatArray ->
                    when {
                        aAndBAreFlipped -> resultPlaneCallback.V[_I, queryBIndex] = it
                        else -> resultPlaneCallback.V[queryBIndex] = it
                    }
                }
            }

        var bestAIndex: Int = -1
        var bestBIndex: Int = -1

        var bestCorrelation = reducer.initializer
        var previousQueryB: F64FlatArray? = null
        TB.mapWindowedWithIndex(m) { queryBIndex, queryB ->

            // skip for continuity
            if (queryB.any { it.isNaN() }) {
                considerCacheQueryBTimeSeriesASlidingDotProductsEmpty.value = true

                return@mapWindowedWithIndex
            }

            var aIndex: Int = -1

            val lineResult = massPearson(
                timeSeries = TANoNaN,
                query = queryB,
                reducer = reducer,

                // vertical speed up
                queryMean = timeSeriesBSlidingMeans[queryBIndex],
                queryStd = timeSeriesBSlidingStds[queryBIndex],
                timeSeriesSlidingMeans = timeSeriesASlidingMeans,
                timeSeriesSlidingStds = timeSeriesASlidingStds,

                // horizontal speed up
                considerCacheQueryTimeSeriesSlidingDotProductsEmpty = considerCacheQueryBTimeSeriesASlidingDotProductsEmpty,
                cacheQueryTimeSeriesSlidingDotProducts = cacheQueryBTimeSeriesASlidingDotProducts,
                previousQuery = previousQueryB,

                resultIndexCallback = when (resultIndexCallback) {
                    null -> null
                    else -> {
                        { aI ->
                            aIndex = aI
                        }
                    }
                },
                resultLineCallback = getResultLineCallback(queryBIndex),
            )

            if (reducer.firstIsBetterThanSecond(first = lineResult, second = bestCorrelation)) {
                bestCorrelation = reducer(bestCorrelation, lineResult)

                @Suppress("SENSELESS_COMPARISON")
                if (resultIndexCallback != null) { // don't listen to this hint, it's wrong
                    bestAIndex = aIndex
                    bestBIndex = queryBIndex
                }
            }
            previousQueryB = queryB

        }

        resultIndexCallback?.invoke(bestAIndex, bestBIndex)

        return bestCorrelation
    }

    /**
     * The MASS-Pearson function ([DoubleArray] variant). It finds the best Pearson correlation value between [query] and all the windows in
     * [timeSeries] which can optionally be aggregated by the queries in [aggregationQueries] at each window position.
     *
     * @param timeSeries the time series given as a whole
     * @param query a query of another time series
     * @param aggregationQueries (optional) should be 1) a query-sized aggregation flat array, or 2) an array of query
     *  sized aggregation flat arrays
     * @param reducer (optional, default = MAX) either [MAX] or [MIN]
     * @param queryDotAggregationProducts (optional) [Double] which can consist of the result of aggregating the
     *  [aggregationQueries] by element wise product and using array that to take the dot product with [query]
     * @param queryTimeSeriesSlidingDotProducts (optional) an array containing the sliding dot products between [query]
     *  and [timeSeries] so that it doesn't have to be calculated in the function anymore.
     * @param queryMean (optional) the mean of [query], so that it doesn't have to be calculated in the function anymore.
     * @param queryStd (optional) the standard deviation of [query], so that it doesn't have to be calculated in the function anymore.
     * @param timeSeriesSlidingMeans (optional) the sliding means of [timeSeries], so that it doesn't have to be
     *  calculated in the function anymore. Don't provide if providing [precalculatedTimeSeriesSlidingMeans].
     * @param timeSeriesSlidingStds (optional) the sliding standard deviations of [timeSeries], so that it doesn't have
     *  to be calculated in the function anymore. Don't provide if providing [precalculatedTimeSeriesSlidingMeans].
     * @param precalculatedTimeSeriesSlidingMeans (optional) the sliding means of [timeSeries], provided to calculate
     *  sliding means and standard deviations with aggregations faster. Don't provide if providing [timeSeriesSlidingMeans]
     *  and [timeSeriesSlidingStds].
     *
     * @param considerCacheQueryTimeSeriesSlidingDotProductsEmpty (optional) delegate boolean to indicate whether [cacheQueryTimeSeriesSlidingDotProducts] should be considered empty.
     * @param cacheQueryTimeSeriesSlidingDotProducts (optional) sliding dot products cache between [query] and [timeSeries]. Needs [considerCacheQueryTimeSeriesSlidingDotProductsEmpty] and [previousQuery].
     * @param previousQuery (optional) the query of the previous iteration, can be [null].
     *
     * @param resultLineCallback (optional) a callback to get all the resulting correlations instead of just the best.
     * @param resultIndexCallback (optional) a callback to get the index of the start of the window position in [timeSeries] for the best found correlation.
     *
     * @return the value of the best Pearson correlation.
     */
    @JvmOverloads
    open fun massPearson(
        timeSeries: DoubleArray,
        query: DoubleArray,
        aggregationQueries: Array<DoubleArray>? = null, // can be supplied to provide aggregation to timeSeries
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        // optional vertical speedup data
        queryDotAggregationProducts: Double? = null, // can be supplied to prevent double calculations
        queryTimeSeriesSlidingDotProducts: DoubleArray? = null,
        queryMean: Double? = null,
        queryStd: Double? = null,
        timeSeriesSlidingMeans: DoubleArray? = null,
        timeSeriesSlidingStds: DoubleArray? = null,
        precalculatedTimeSeriesSlidingMeans: DoubleArray? = null, // for calculating sliding mean/std with aggregation faster

        // optional horizontal speedup data
        considerCacheQueryTimeSeriesSlidingDotProductsEmpty: Delegate<Boolean>? = null,
        cacheQueryTimeSeriesSlidingDotProducts: DoubleArray? = null,
        previousQuery: DoubleArray? = null,

        // get optional result line
        resultLineCallback: ((DoubleArray) -> Unit)? = null,

        // get optional result index
        resultIndexCallback: ((Int) -> Unit)? = null,
    ): Double = massPearson(
        timeSeries = timeSeries.asF64Array(),
        query = query.asF64Array(),
        aggregationQueries = aggregationQueries?.let {
            F64Array(it.size, it.first().size) { i, j -> it[i][j] }
        },
        reducer = reducer,
        queryDotAggregationProducts = queryDotAggregationProducts,
        queryTimeSeriesSlidingDotProducts = queryTimeSeriesSlidingDotProducts?.asF64Array(),
        queryMean = queryMean,
        queryStd = queryStd,
        timeSeriesSlidingMeans = timeSeriesSlidingMeans?.asF64Array(),
        timeSeriesSlidingStds = timeSeriesSlidingStds?.asF64Array(),
        precalculatedTimeSeriesSlidingMeans = precalculatedTimeSeriesSlidingMeans?.asF64Array(),
        considerCacheQueryTimeSeriesSlidingDotProductsEmpty = considerCacheQueryTimeSeriesSlidingDotProductsEmpty,
        cacheQueryTimeSeriesSlidingDotProducts = cacheQueryTimeSeriesSlidingDotProducts?.asF64Array(),
        previousQuery = previousQuery?.asF64Array(),
        resultLineCallback = resultLineCallback?.let {
            {
                resultLineCallback(it.toDoubleArray())
            }
        },
        resultIndexCallback = resultIndexCallback,
    )

    /**
     * The MASS-Pearson function. It finds the best Pearson correlation value between [query] and all the windows in
     * [timeSeries] which can optionally be aggregated by the queries in [aggregationQueries] at each window position.
     *
     * @param timeSeries the time series given as a whole
     * @param query a query of another time series
     * @param aggregationQueries (optional) should be 1) a query-sized aggregation flat array, or 2) an array of query
     *  sized aggregation flat arrays
     * @param reducer (optional, default = MAX) either [MAX] or [MIN]
     * @param queryDotAggregationProducts (optional) [Double] which can consist of the result of aggregating the
     *  [aggregationQueries] by element wise product and using array that to take the dot product with [query]
     * @param queryTimeSeriesSlidingDotProducts (optional) an array containing the sliding dot products between [query]
     *  and [timeSeries] so that it doesn't have to be calculated in the function anymore.
     * @param queryMean (optional) the mean of [query], so that it doesn't have to be calculated in the function anymore.
     * @param queryStd (optional) the standard deviation of [query], so that it doesn't have to be calculated in the function anymore.
     * @param timeSeriesSlidingMeans (optional) the sliding means of [timeSeries], so that it doesn't have to be
     *  calculated in the function anymore. Don't provide if providing [precalculatedTimeSeriesSlidingMeans].
     * @param timeSeriesSlidingStds (optional) the sliding standard deviations of [timeSeries], so that it doesn't have
     *  to be calculated in the function anymore. Don't provide if providing [precalculatedTimeSeriesSlidingMeans].
     * @param precalculatedTimeSeriesSlidingMeans (optional) the sliding means of [timeSeries], provided to calculate
     *  sliding means and standard deviations with aggregations faster. Don't provide if providing [timeSeriesSlidingMeans]
     *  and [timeSeriesSlidingStds].
     *
     * @param considerCacheQueryTimeSeriesSlidingDotProductsEmpty (optional) delegate boolean to indicate whether [cacheQueryTimeSeriesSlidingDotProducts] should be considered empty.
     * @param cacheQueryTimeSeriesSlidingDotProducts (optional) sliding dot products cache between [query] and [timeSeries]. Needs [considerCacheQueryTimeSeriesSlidingDotProductsEmpty] and [previousQuery].
     * @param previousQuery (optional) the query of the previous iteration, can be [null].
     *
     * @param resultLineCallback (optional) a callback to get all the resulting correlations instead of just the best.
     * @param resultIndexCallback (optional) a callback to get the index of the start of the window position in [timeSeries] for the best found correlation.
     *
     * @return the value of the best Pearson correlation.
     */
    @JvmOverloads
    open fun massPearson(
        timeSeries: F64FlatArray,
        query: F64FlatArray,
        aggregationQueries: F64Array? = null, // can be supplied to provide aggregation to timeSeries
        reducer: AggregationWithReducerWithArgWithCompare = MAX,

        // optional vertical speedup data
        queryDotAggregationProducts: Double? = null, // can be supplied to prevent double calculations
        queryTimeSeriesSlidingDotProducts: F64FlatArray? = null,
        queryMean: Double? = null,
        queryStd: Double? = null,
        timeSeriesSlidingMeans: F64FlatArray? = null,
        timeSeriesSlidingStds: F64FlatArray? = null,
        precalculatedTimeSeriesSlidingMeans: F64FlatArray? = null, // for calculating sliding mean/std with aggregation faster

        // optional horizontal speedup data
        considerCacheQueryTimeSeriesSlidingDotProductsEmpty: Delegate<Boolean>? = null,
        cacheQueryTimeSeriesSlidingDotProducts: F64FlatArray? = null,
        previousQuery: F64FlatArray? = null,

        // get optional result line
        resultLineCallback: ((F64FlatArray) -> Unit)? = null,

        // get optional result index
        resultIndexCallback: ((Int) -> Unit)? = null,
    ): Double {
        val T = timeSeries
        val Q = query
        val QTCache = cacheQueryTimeSeriesSlidingDotProducts

        val n = T.length
        val m = Q.length

        require(T.all { it.isNotNaN() }) { "timeSeries cannot contain NaNs for removed continuity" }

        // should be an array of arrays of size m or one array of size m
        val Qa = aggregationQueries?.let {
            if (it.nDim == 1) {
                require(it.length == m) { "aggregation should be same size as query" }
                return@let it.reshape(1, m)
            } else {
                require(it.nDim == 2) { "aggregation-array should only consist of flat arrays" }
                require(it.shape[1] == m) { "aggregation-array should have flat arrays of query size" }
                return@let it
            }
        }

        // number of queries (including aggregation queries)
        val q = (Qa?.length ?: 0) + 1.0

        val QDotsQa = queryDotAggregationProducts ?: Qa?.let { _ ->
            arrayCache[Qa.shape[1], Slot.R].also {
                for (i in it.indices)
                    it[i] = Qa.view(i, 1).asFlatArray().product()
            } dot Q
        }

        val QT = when {
            queryTimeSeriesSlidingDotProducts != null ->
                queryTimeSeriesSlidingDotProducts

            QTCache != null && considerCacheQueryTimeSeriesSlidingDotProductsEmpty != null ->
                QTCache.also {
                    updateQTCache(
                        QTCache = it,
                        considerQTCacheEmpty = considerCacheQueryTimeSeriesSlidingDotProductsEmpty,
                        Q = Q,
                        T = T,
                        n = n,
                        m = m,
                        previousQuery = previousQuery,
                        returnArraySlot = Slot.B3,
                    )
                }

            else ->
                slidingDotProducts(Q = Q, T = T, returnArraySlot = Slot.B3)
        }

        // apply aggregation for T on QT
        if (Qa != null && QDotsQa != null) {
            QT += QDotsQa
            QT /= q
        }

        // compute mean and stds
        val mu_Q = (queryMean ?: Q.mean())
            .let { if (it.isInfinite() || it.isNaN()) 0.0 else it }

        val sigma_Q = (queryStd ?: Q.let {
            it -= mu_Q
            sqrt(it.dot() / it.length)
        })
            .let { if (it.isInfinite() || it.isNaN()) 0.0 else it }

        val mus_T: F64FlatArray
        val sigmas_T: F64FlatArray
        if (timeSeriesSlidingMeans != null && timeSeriesSlidingStds != null) {
            mus_T = arrayCache[timeSeriesSlidingMeans.length, Slot.O]
            timeSeriesSlidingMeans.copyTo(mus_T)

            sigmas_T = timeSeriesSlidingStds
        } else {
            mus_T = arrayCache[n - m + 1, Slot.O]
            sigmas_T = arrayCache[n - m + 1, Slot.P]

            if (Qa == null) {
                computeSlidingMeanStd(
                    T = T,
                    windowSize = m,
                    destMeanArray = mus_T,
                    destStdArray = sigmas_T,
                )
            } else {
                computeSlidingMeanStdWithAgg(
                    T = T,
                    windowSize = m,
                    destMeanArray = mus_T,
                    destStdArray = sigmas_T,
                    aggregations = Qa,
                    precalculatedTMeans = precalculatedTimeSeriesSlidingMeans,
                )
            }
        }

        mus_T.transformInPlace { if (it.isInfinite() || it.isNaN()) 0.0 else it }
        sigmas_T.transformInPlace { if (it.isInfinite() || it.isNaN()) 0.0 else it }

        // calculate the pearson correlations
        // formula:
        //       (QT - (m * mu_Q * mus_T)) /
        //            (m * sigma_Q * sigmas_T)

        val P = calculatePearsonProfile(
            m = m,
            QT = QT,
            mu_Q = mu_Q,
            mus_T = mus_T,
            sigma_Q = sigma_Q,
            sigmas_T = sigmas_T
        )

        // fix potential divide by 0 issues
        P.transformInPlace { if (it.isInfinite() || it.isNaN()) 0.0 else it }

        // reverse aggregation for T on QT
        if (Qa != null && QDotsQa != null) {
            QT *= q
            QT -= QDotsQa
        }

        resultLineCallback?.invoke(P)

        val index = reducer.argFunction(P)
        resultIndexCallback?.invoke(index)

        if (P.any { it !in -1.01..1.01 })
            throw IllegalArgumentException("$P")

        return P[index]
    }

    /**
     * Updates the given [QTCache] for the next horizontal iteration.
     *
     * @param QTCache the cache to update.
     * @param considerQTCacheEmpty a delegate to a boolean value to indicate whether [QTCache] should be considered
     * empty and over-writable.
     * @param Q the current query.
     * @param T the current time series.
     * @param n (optional) the size of [T].
     * @param m (optional) the size of [Q] (or the window size).
     * @param previousQuery the previous query, can be [null] if at the first iteration.
     * @param returnArraySlot (optional) the array from the cache that is used in the [slidingDotProducts] call if the cache is empty.
     */
    @Suppress("UNUSED_VALUE")
    @JvmOverloads
    open fun updateQTCache(
        QTCache: F64FlatArray,
        considerQTCacheEmpty: Delegate<Boolean>,
        Q: F64FlatArray,
        T: F64FlatArray,
        n: Int = T.length,
        m: Int = Q.length,
        previousQuery: F64FlatArray?,
        returnArraySlot: Slot = Slot.B,
    ) {
        var cacheIsEmpty by considerQTCacheEmpty

        if (cacheIsEmpty) {
            slidingDotProducts(Q = Q, T = T, returnArraySlot = returnArraySlot).copyTo(QTCache)
            cacheIsEmpty = false
        } else {
            val prevQFirst = previousQuery!![0]
            val Qlast = Q[m - 1]

            for (i in n - m downTo 1) {
                QTCache[i] = QTCache[i - 1] - prevQFirst * T[i - 1] + Qlast * T[i + m - 1]
            }
            QTCache[0] = Q dot T.sliceFlat(to = m)
        }
    }

    /**
     * Calculation of sliding dot products.
     * @param Q a query
     * @param T a time series
     * @param n (optional) the length of [T]
     * @param m (optional) the length of query [Q] (or the windowSize)
     * @param returnArraySlot (optional) the array slot in the cache that will be used as return value
     * @return the sliding dot products between [Q] and [T]
     */
    @JvmOverloads
    fun slidingDotProducts(
        Q: F64FlatArray,
        T: F64FlatArray,
        n: Int = T.length,
        m: Int = Q.length,
        returnArraySlot: Slot = Slot.B,
    ): F64FlatArray {
        // using array cache to save creating new arrays all the time.
        val nAppended = nextPowerOf2(n)

        val TAppended = arrayCache[2 * nAppended, Slot.A] { 0.0 }
        T.copyTo(TAppended.sliceFlat(from = nAppended - n, to = nAppended))

        val fourier: DoubleFFT_1D = fourierCache.getOrPut(nAppended) { DoubleFFT_1D(nAppended.toLong()) }

        // reverse Q to Qr and append Qr with zeroes and double size for in-place fft
        val Qra = arrayCache[2 * nAppended, returnArraySlot] {
            if (it < m) Q[m - it - 1] else 0.0
        }

        // perform fft on Q_ra
        val Qraf = Qra.also { fourier.realForwardFull(it.data) }

        val Taf = TAppended.also { fourier.realForwardFull(it.data) }

        // !! overwrites [Q_raf]
        val multiplication = Qraf.also {
            for (k in 0 until it.length / 2) {
                val aRe = it[2 * k]
                val aIm = it[2 * k + 1]
                val bRe = Taf[2 * k]
                val bIm = Taf[2 * k + 1]

                it[2 * k] = aRe * bRe - aIm * bIm
                it[2 * k + 1] = aRe * bIm + aIm * bRe
            }
        }

        val QT2 = multiplication.also { fourier.complexInverse(it.data, true) }

        return QT2
            .sliceFlat(step = 2) // take only reals
            .sliceFlat(from = m - 1 + (nAppended - n)) // remove padding
    }

    /**
     * Efficient moving average and moving std calculations with support of any number of aggregation queries.
     * result is stored in destMeanArray and destStdArray (which must be size = T.size - windowSize + 1)
     * Compatible with continuity removed (NaNs in [T]), but it is better to give it the original time series (without NaNs) for speed.
     *
     * @param T Time series over which to slide
     * @param windowSize the window size
     * @param aggregations array of aggregation queries.
     * @param precalculatedTMeans Optional speed-up
     *
     * @return sliding means and -standard deviations with aggregations
     */
    @JvmOverloads
    open fun computeSlidingMeanStdWithAgg(
        T: F64FlatArray,
        windowSize: Int,
        aggregations: F64Array,
        precalculatedTMeans: F64FlatArray? = null,
    ): Pair<F64FlatArray, F64FlatArray> {
        val size = T.length
        val meanArray = F64FlatArray(size - windowSize + 1)
        val stdArray = F64FlatArray(size - windowSize + 1)

        computeSlidingMeanStdWithAgg(
            T = T,
            windowSize = windowSize,
            aggregations = aggregations,
            destMeanArray = meanArray,
            destStdArray = stdArray,
            precalculatedTMeans = precalculatedTMeans,
        )

        return meanArray to stdArray
    }

    /**
     * Efficient moving average and moving std calculations with support of any number of aggregation queries.
     * result is stored in destMeanArray and destStdArray (which must be size = T.size - windowSize + 1)
     * Compatible with continuity removed (NaNs in [T]), but it is better to give it the original time series (without NaNs) for speed.
     *
     * @param T Time series over which to slide
     * @param windowSize the window size
     * @param aggregations array of aggregation queries.
     * @param destMeanArray the destination array for the means, must be `T.size - windowSize + 1`
     * @param destStdArray the destination array for the stds, must be `T.size - windowSize + 1`
     * @param precalculatedTMeans Optional speed-up
     */
    @JvmOverloads
    open fun computeSlidingMeanStdWithAgg(
        T: F64FlatArray,
        windowSize: Int,
        aggregations: F64Array,
        destMeanArray: F64FlatArray,
        destStdArray: F64FlatArray,
        precalculatedTMeans: F64FlatArray? = null,
    ) {
        if (T.any { it.isNaN() }) println("WARNING: slidingMeanAndStd is called with a continuity-removed time series. Results might not be as fast.")
        val n = T.length
        val m = windowSize
        require(destMeanArray.length == n - m + 1)
        require(destStdArray.length == n - m + 1)
        var Qa = aggregations

        if (Qa.nDim == 1) {
            require(Qa.length == m) { "aggregation should be window size" }
            Qa = Qa.reshape(1, m)
        } else {
            require(Qa.nDim == 2) { "aggregation-array should only consist of flat arrays" }
            require(Qa.shape[1] == m) { "aggregation-array should have flat arrays of query size" }
        }

        val q = Qa.length.toDouble() + 1.0
        val qSquared = q * q

        // sum of means of Qas
        var QaMeansSum = 0.0
        for (it in 0 until Qa.length) {
            QaMeansSum += Qa.viewAsFlat(it).mean()
        }

        // Qas aggregated by sum
        val QaSums = arrayCache[m, Slot.Q] {
            Qa.viewAsFlat(it, 1).sum()
        }

        val QaSumsSum = QaSums.sum()

        val dotOfQaSums = QaSums.dot()

        // first calculate means using buffer method since aggregation can be done afterwards
        var sumBuffer = Double.NaN
        var squareBuffer = Double.NaN

        val TQaSumsSlidingDots =
            slidingDotProducts(Q = QaSums, T = T, returnArraySlot = Slot.B4)

        if (precalculatedTMeans != null) {
            precalculatedTMeans.copyTo(destMeanArray)
            destMeanArray += QaMeansSum
            destMeanArray /= q
        }

        for (i in 0..n - m) {
            if (sumBuffer.isNaN() || squareBuffer.isNaN()) {
                T.sliceFlat(from = i, to = i + m).run {
                    sumBuffer = sum()
                    squareBuffer = dot()
                }
            }
            if (sumBuffer.isNaN() || squareBuffer.isNaN()) {
                destMeanArray[i] = Double.NaN
                destStdArray[i] = Double.NaN
                continue
            }

            var aggMean: Double
            if (precalculatedTMeans == null) {
                val mean = sumBuffer / m
                aggMean = (mean + QaMeansSum) / q
                destMeanArray[i] = aggMean
            } else {
                aggMean = destMeanArray[i]
            }

            var stdSum = 0.0

            stdSum += squareBuffer / qSquared

            stdSum += 2.0 * TQaSumsSlidingDots[i] / qSquared // stdSum += 2.0 * (Tw * QaSums).sum() / qSquared

            stdSum -= 2.0 * aggMean * sumBuffer / q // stdSum -= 2.0 * aggMean * Tw.sum() / q

            stdSum += dotOfQaSums / qSquared

            stdSum -= 2.0 * aggMean * QaSumsSum / q // stdSum -= 2.0 * aggMean * QaSums.sum() / q

            stdSum += m * (aggMean * aggMean)

            destStdArray[i] = sqrt(stdSum / m.toDouble())

            if (i < n - m) {
                val toBeRemoved = T[i]
                val toBeAdded = T[i + m]

                sumBuffer += toBeAdded - toBeRemoved
                squareBuffer += toBeAdded * toBeAdded - toBeRemoved * toBeRemoved
            }
        }
    }

    // These functions don't need any cache, so they are made statically available using Kotlin's companion object.
    // Access from Java using `StampPearson.Companion.function()`
    // Access from Kotlin using `StampPearson.function()`
    companion object {

        /**
         * Can convert an array containing Euclidean distances to Pearson correlations.
         */
        fun euclideanToPearsonCorrelations(series: F64FlatArray, windowSize: Int): F64FlatArray =
            (1.0 - (1.0 / (2.0 * windowSize)) * (series * series)).asFlatArray()

        /**
         * Can convert an array containing Pearson correlations to Euclidean distances.
         */
        fun pearsonCorrelationsToEuclidean(series: F64FlatArray, windowSize: Int): F64FlatArray =
            (2.0 * windowSize * (1.0 - series)).asFlatArray()

        /**
         * Extremely efficient moving average and moving std calculations.
         * result is stored in destMeanArray and destStdArray (which must be size = T.size - windowSize + 1)
         * Compatible with continuity removed (NaNs in [T]), but it is better to give it the original time series (without NaNs) for speed.
         * @see [computeSlidingMeanStd]
         *
         * T. Rakthanmanon et al., “Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping,”
         *
         * @param T Time series over which to slide
         * @param windowSize the window size
         *
         * @return sliding means and -standard deviations
         */
        fun computeSlidingMeanStd(T: F64FlatArray, windowSize: Int): Pair<F64FlatArray, F64FlatArray> {
            val size = T.length
            val meanArray = F64FlatArray(size - windowSize + 1)
            val stdArray = F64FlatArray(size - windowSize + 1)
            computeSlidingMeanStd(T, windowSize, meanArray, stdArray)

            return meanArray to stdArray
        }

        /**
         * Extremely efficient moving average and moving std calculations.
         * result is stored in destMeanArray and destStdArray (which must be size = T.size - windowSize + 1)
         * Compatible with continuity removed (NaNs in [T]), but it is better to give it the original time series (without NaNs) for speed.
         * @see [computeSlidingMeanStd]
         *
         * T. Rakthanmanon et al., “Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping,”
         *
         * @param T Time series over which to slide
         * @param windowSize the window size
         * @param destMeanArray the destination array for the means, must be `T.size - windowSize + 1`
         * @param destStdArray the destination array for the stds, must be `T.size - windowSize + 1`
         */
        fun computeSlidingMeanStd(
            T: F64FlatArray,
            windowSize: Int,
            destMeanArray: F64FlatArray,
            destStdArray: F64FlatArray,
        ) {
            if (T.any { it.isNaN() }) println("WARNING: slidingMeanAndStd is called with a continuity-removed time series. Results might not be as fast.")
            val n = T.length
            val m = windowSize
            require(destMeanArray.length == n - m + 1)
            require(destStdArray.length == n - m + 1)

            var sumBuffer = Double.NaN
            var squareBuffer = Double.NaN

            for (i in 0..n - m) {
                if (sumBuffer.isNaN() || squareBuffer.isNaN()) {
                    T.sliceFlat(from = i, to = i + m).run {
                        sumBuffer = sum()
                        squareBuffer = dot()
                    }
                }
                if (sumBuffer.isNaN() || squareBuffer.isNaN()) {
                    destMeanArray[i] = Double.NaN
                    destStdArray[i] = Double.NaN
                    continue
                }

                val mean = sumBuffer / m
                destMeanArray[i] = mean

                val square = mean * mean
                destStdArray[i] = sqrt(squareBuffer / m - square)

                if (i < n - m) {
                    val toBeRemoved = T[i]
                    val toBeAdded = T[i + m]

                    sumBuffer += toBeAdded - toBeRemoved
                    squareBuffer += toBeAdded * toBeAdded - toBeRemoved * toBeRemoved
                }
            }
        }

        /**
         * Calculate Pearson efficiently using matrix multiplications.
         * Formula:
         *       (QT - (m * mu_Q * mus_T)) /
         *            (m * sigma_Q * sigmas_T)
         * @param m window size
         * @param QT dot products between query Q and all possible positions in T (will be overwritten!!)
         * @param mu_Q mean of Q
         * @param mus_T means of all possible placements of a window in T (will be overwritten!!)
         * @param sigma_Q standard deviation of Q
         * @param sigmas_T standard deviation of all possible placements of a window in T
         * @return sliding Pearson correlation measure array.
         */
        fun calculatePearsonProfile(
            m: Int,
            QT: F64FlatArray, // will not be overwritten
            mu_Q: Double,
            mus_T: F64FlatArray, // will be overwritten
            sigma_Q: Double,
            sigmas_T: F64FlatArray, // will not be overwritten
        ): F64FlatArray {
            // !! overwriting [mus_T]
            val partNominator = mus_T.also { it *= mu_Q * -m }

            // !! overwriting [partNominator]
            val nominator = partNominator.also { it += QT }

            // !! overwriting [nominator]
            val answerWithoutSigmasT = nominator.also { it /= m.toDouble() * sigma_Q }

            // !! overwriting [answerWithoutSigmasT]
            val answer = answerWithoutSigmasT.also { it /= sigmas_T }

            return answer
        }
    }
}
