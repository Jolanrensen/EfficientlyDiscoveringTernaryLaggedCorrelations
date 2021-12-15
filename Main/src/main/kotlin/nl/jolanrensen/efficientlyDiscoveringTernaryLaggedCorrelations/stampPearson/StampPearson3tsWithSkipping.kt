@file:Suppress("LocalVariableName")

package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.AggregationWithReducerWithArgWithCompare
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.Continuity.numberOfWindowPlacementsSaved
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3tsWithSkipping.AggregationMethod.NONE
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3tsWithSkipping.TimeSeries.*
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampArrayCache.Slot
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor._I
import org.jtransforms.fft.DoubleFFT_1D
import java.io.Serializable
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * TODO
 */
class StampPearson3tsWithSkipping @JvmOverloads constructor(
    val continuitySkipping: ContinuitySkippingType,
    arrayCache: StampArrayCache = StampArrayCache(),
    fourierCache: MutableMap<Int, DoubleFFT_1D> = mutableMapOf(),
) : StampPearson(arrayCache, fourierCache), Serializable {

    @JvmOverloads
    constructor(
        maxArraySize: Int,
        continuitySkipping: ContinuitySkippingType,
        slotManager: StampArrayCache.SlotManager = StampArrayCache.SlotsWithReuse,
    ) : this(
        continuitySkipping = continuitySkipping,
        arrayCache = StampArrayCache(
            maxSize = maxArraySize,
            slotManager = slotManager,
        ),
    )

    companion object {

        /**
         * @param
         * @param endOfSliceIndex exclusive index of end of slice of time series
         * @return function that can be executed to reverse the fix applied to [timeSeries], `null` if nothing changed.
         */
        fun fixRemovedContinuityWithLagBound(
            timeSeries: F64FlatArray,
            timeSeriesNoNaN: F64FlatArray,
            endOfSliceIndex: Int,
            windowSize: Int,
            continuitySkipping: ContinuitySkippingType,
        ): (() -> Unit)? { // TODO this only works for straight lines, not MP

            if (continuitySkipping == ContinuitySkippingType.MATRIX_PROFILE)
                throw IllegalArgumentException("unsupported")

            if (endOfSliceIndex >= timeSeries.length) {
                // no fixing necessary
                return null
            }

            fun addValueToLeftSideOfNaNs(indexInNaNs: Int): () -> Unit {
                // find left of NaNs
                var j = indexInNaNs
                while (timeSeries[--j].isNaN() && j > 0) Unit
                j++

//            // find right of NaNs (exclusive)
//            var k = indexInNaNs
//            while (timeSeries[++k].isNaN() && k < timeSeries.size - 1) Unit
//
//            // add one more value on left side of NaNs based on the slope on the right side
//            val slope = timeSeries[k + 1] - timeSeries[k]
//            timeSeries[j] = (k - j + 1) * slope
                timeSeries[j] = timeSeriesNoNaN[j]

                return {
                    timeSeries[j] = Double.NaN
                }
            }

            if (timeSeries[endOfSliceIndex - 1].isNaN()) {
                if (timeSeries.sliceFlat(from = endOfSliceIndex - 1, to = timeSeries.length).all { it.isNaN() }) {
                    // no fixing necessary
                    return null
                }

                return addValueToLeftSideOfNaNs(endOfSliceIndex - 1)
            }

            // if end of slice is next to NaN-piece in the middle of a time series
            for ((i, value) in timeSeries.sliceFlat(
                from = (endOfSliceIndex - windowSize).coerceAtLeast(0),
                to = endOfSliceIndex).withIndex()
            ) {
                if (value.isNaN()) {
                    return addValueToLeftSideOfNaNs(i + 1)
                }
            }

            return null
        }
    }

    enum class ContinuitySkippingType {
        MATRIX_PROFILE, STRAIGHT_LINE
    }

    enum class TimeSeries {
        A, B, C
    }

    enum class AggregationMethod(val textual: String, val aggregatedTimeSeries: List<TimeSeries>) :
        Serializable {
        AB("agg(a, b) and c", listOf(A, B)),
        BC("a and agg(b, c)", listOf(B, C)),
        AC("agg(a, c) and b", listOf(A, C)),
        NONE("Aggregation didn't happen", emptyList()),
    }

    data class SkippableIndex(
        val aggregationMethod: AggregationMethod,
        val aIndex: Int,
        val bIndex: Int,
        val cIndex: Int,
    )

    /**
     * Find the best correlation between 3 time series for each possible window position in each.
     *
     * @param timeSeriesA First time series (must have the smallest amount of `NaN`s if any).
     * @param timeSeriesB Second time series (must have the second to smallest amount of `NaN`s if any).
     * @param timeSeriesC Third time series (will be split first) (must have the largest amount of `NaN`s if any).
     * @param windowSize, also called `m`, the window size used to traverse all three time series.
     * @param reducer Defines whether to find the highest or lowest correlation.
     *      In practice this is either [MIN] or [MAX]. Defaults to [MAX].
     *
     * @param timeSeriesAWithoutNaN Optional, but if any continuity is removed from [timeSeriesA], this needs to be supplied.
     * @param timeSeriesBWithoutNaN Optional, but if any continuity is removed from [timeSeriesB], this needs to be supplied.
     * @param timeSeriesCWithoutNaN Optional, but if any continuity is removed from [timeSeriesC], this needs to be supplied.
     *
     * @param timeSeriesAWithoutNaNSlidingMeans Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
     *      [timeSeriesA]. NOTE: As the name suggests, if [timeSeriesA] contains `NaN`s,
     *      use [timeSeriesAWithoutNaN] to calculate it.
     * @param timeSeriesAWithoutNaNSlidingStds Optional optimization. Can be calculated using [StampPearson.computeSlidingMeanStd] of
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
     * @param overwrite If `true`, then the provided time series might be written to, else they are copied. Defaults to `false`.
     * @param useQTCache If `true`, then the horizontal sliding dots cache optimization is used. Defaults to `true`.
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

        timeSeriesAWithoutNaN: F64FlatArray? = null,
        timeSeriesBWithoutNaN: F64FlatArray? = null,
        timeSeriesCWithoutNaN: F64FlatArray? = null,

        timeSeriesAWithoutNaNSlidingMeans: F64FlatArray? = null,
        timeSeriesAWithoutNaNSlidingStds: F64FlatArray? = null,
        timeSeriesBSlidingMeans: F64FlatArray? = null,
        timeSeriesBSlidingStds: F64FlatArray? = null,
        timeSeriesCSlidingMeans: F64FlatArray? = null,
        timeSeriesCSlidingStds: F64FlatArray? = null,

//        overwrite: Boolean = false,
//        useQTCache: Boolean = true,
        lagBound: Int? = null,
//        skipIndices: List<SkippableIndex>? = null,

        // get optional result matrix
        resultCubeCallback: F64Array? = null,
        resultIndexCallback: ((aggregation: AggregationMethod, aIndex: Int, bIndex: Int, cIndex: Int) -> Unit)? = null,
    ): Double {
        val m = windowSize

        val TA: F64FlatArray = timeSeriesA
        val TB: F64FlatArray = timeSeriesB
        val TC: F64FlatArray = timeSeriesC


//        var timeSeriesAWithoutNaN = timeSeriesAWithoutNaN
//        var timeSeriesBWithoutNaN = timeSeriesBWithoutNaN
//        var timeSeriesCWithoutNaN = timeSeriesCWithoutNaN


        val nanCountA = numberOfWindowPlacementsSaved(TA, m)
        val nanCountB = numberOfWindowPlacementsSaved(TB, m)
        val nanCountC = numberOfWindowPlacementsSaved(TC, m)

        if (nanCountA > 0) require(timeSeriesAWithoutNaN != null) { "The supplied time series A has continuity removed, please supply the original time series without NaNs as well." }
        if (nanCountB > 0) require(timeSeriesBWithoutNaN != null) { "The supplied time series B has continuity removed, please supply the original time series without NaNs as well." }
        if (nanCountC > 0) require(timeSeriesCWithoutNaN != null) { "The supplied time series C has continuity removed, please supply the original time series without NaNs as well." }

        require(nanCountA <= nanCountB && nanCountB <= nanCountC) { "timeSeriesA must have the least NaNs, then timeSeriesB, then timeSeriesC" }

        val TANoNaN = timeSeriesAWithoutNaN ?: TA
        val TBNoNaN = timeSeriesBWithoutNaN ?: TB
        val TCNoNaN = timeSeriesCWithoutNaN ?: TC

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

        var timeSeriesASlidingMeans = timeSeriesAWithoutNaNSlidingMeans
        var timeSeriesASlidingStds = timeSeriesAWithoutNaNSlidingStds

        if (timeSeriesASlidingMeans == null || timeSeriesASlidingStds == null) {
            timeSeriesASlidingMeans = arrayCache[nA - m + 1, Slot.F]
            timeSeriesASlidingStds = arrayCache[nA - m + 1, Slot.G]
            computeSlidingMeanStd(
                T = TANoNaN,
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
                T = TBNoNaN,
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
                T = TCNoNaN,
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
                timeSeriesA = TANoNaN,
                timeSeriesB = TB,
                timeSeriesBWithoutNaN = TBNoNaN,
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
        timeSeriesAWithoutNaN: F64FlatArray? = null, // if any continuity is removed from [timeSeriesA], this needs to be supplied
        timeSeriesBWithoutNaN: F64FlatArray? = null, // if any continuity is removed from [timeSeriesB], this needs to be supplied
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

//        require(TA.all { it.isNotNaN() }) { "timeSeriesA is not allowed to have continuity removed and thus may not contain NaNs." }

        val m = queryC.length
        val nA = TA.length
        var nB = TB.length

        if (TA.any { it.isNaN() }) require(timeSeriesAWithoutNaN != null) { "The supplied time series A has continuity removed, please supply the original time series without NaNs as well." }
        if (TB.any { it.isNaN() }) require(timeSeriesBWithoutNaN != null) { "The supplied time series B has continuity removed, please supply the original time series without NaNs as well." }

        val TANoNaN = timeSeriesAWithoutNaN ?: TA
        var TBNoNaN = timeSeriesBWithoutNaN ?: TB

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
        var restoreTBContinuityLagBoundFix: (() -> Unit)? = null

        if (lagBound != null) {

            // the bound would make that no m values are left in TB
            if (queryCIndex!! > nB - m) return reducer.initializer

            val startBoundIndex = (queryCIndex - lagBound).coerceAtLeast(0)
            val endBoundIndex = (queryCIndex + lagBound).coerceAtMost(nB - m)

            restoreTBContinuityLagBoundFix = fixRemovedContinuityWithLagBound(
                timeSeries = TB,
                timeSeriesNoNaN = TBNoNaN,
                endOfSliceIndex = endBoundIndex + m,
                windowSize = m,
                continuitySkipping = continuitySkipping,
            )
            TB = TB.sliceFlat(from = startBoundIndex, to = endBoundIndex + m)
            if (TB !== TBNoNaN) TBNoNaN = TBNoNaN.sliceFlat(from = startBoundIndex, to = endBoundIndex + m)


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
                        T = TANoNaN,
                        n = nA,
                        m = m,
                        previousQuery = previousQueryC,
                        returnArraySlot = Slot.B1,
                    )
                }
            else -> slidingDotProducts(
                Q = queryC,
                T = TANoNaN,
                returnArraySlot = Slot.B1,
            )
        }

        // precalculate all QB dot QC
        val queryCTimeSeriesBSlidingDotProducts = when {
            cacheQueryCTimeSeriesBSlidingDotProducts != null && considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty != null ->
                cacheQueryCTimeSeriesBSlidingDotProducts.also {
                    updateQTCache(
                        QTCache = it,
                        considerQTCacheEmpty = considerCacheQueryCTimeSeriesBSlidingDotProductsEmpty,
                        Q = queryC,
                        T = TBNoNaN,
                        n = nB,
                        m = m,
                        previousQuery = previousQueryC,
                        returnArraySlot = Slot.B2,
                    )
                }
            else -> slidingDotProducts(
                Q = queryC,
                T = TBNoNaN,
                returnArraySlot = Slot.B2,
            )
        }

        // reused later in agg(T, Q1) and Q2
        //      where Q1 is from [timeSeriesToBeSplit] and Q2 is [query]
        val queryCMean = queryCMean ?: queryC.mean()
        val queryCStd = queryCStd ?: run {
            arrayCache[queryC.length, Slot.M] { queryC[it] }
                .let {
                    it -= queryCMean
                    sqrt(it.dot() / it.length)
                }
        }

        // reused later in T and agg(Q1, Q2)
        //      where Q1 is from [timeSeriesToBeSplit] and Q2 is [query]
        var timeSeriesASlidingMeans = timeSeriesASlidingMeans
        var timeSeriesASlidingStds = timeSeriesASlidingStds

        if (timeSeriesASlidingStds == null || timeSeriesASlidingMeans == null) {
            timeSeriesASlidingMeans = arrayCache[nA - m + 1, Slot.F]
            timeSeriesASlidingStds = arrayCache[nA - m + 1, Slot.G]
            computeSlidingMeanStd(
                T = TANoNaN, // no difference with TA in practice
                windowSize = m,
                destMeanArray = timeSeriesASlidingMeans,
                destStdArray = timeSeriesASlidingStds,
            )
        }

        val timeSeriesASlidingMeansAggQC = arrayCache[nA - m + 1, Slot.H]
        val timeSeriesASlidingStdsAggQC = arrayCache[nA - m + 1, Slot.I]
        computeSlidingMeanStdWithAgg(
            T = TANoNaN, // no difference with TA in practice
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
                T = TBNoNaN,
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
                timeSeriesA = TANoNaN,
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

        restoreTBContinuityLagBoundFix?.invoke()
        resultIndexCallback?.invoke(bestAggregation, bestAIndex, bestBIndex + bIndexModifier)

        return bestCorrelation
    }

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

        val resultsMatrix = if (null == resultLineCallback) null else F64Array(3, nA - m + 1)
        val resultsIndexArray = if (resultIndexCallback == null) null else IntArray(3) { -1 }



        // agg(TA, QB) and QC
        val a = massPearson(
            timeSeries = TA,
            aggregationQueries = QB,
            query = QC,
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
            query = QB.copy().also {
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
//            precalculatedTimeSeriesSlidingMeans = timeSeriesASlidingMeans, not needed

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
            val aIndex = resultsIndexArray!![aggregation]

            resultIndexCallback(AggregationMethod.values()[aggregation], aIndex + aIndexModifier)
        }


        if (resultLineCallback != null) {
            resultsMatrix!!.transformInPlace { it.coerceNotNaN(reducer.initializer) }
            var i = 0
            resultsMatrix.V[0].transformInPlace { reducer(resultsMatrix.V[_I, i++]) }

            val line = resultsMatrix.V[0].asFlatArray()
            line.copyTo(resultLineCallback.sliceFlat(from = aIndexModifier, to = aIndexModifier + line.length))
        }

        return bestCorrelation
    }


    /**
     * @param aggregationQueries should be 1) a query-sized aggregation flat array, or 2) an array of query sized aggregation flat arrays
     * @param queryDotAggregationProducts is an optional [Double] which can consist of the result of aggregating the
     *  [aggregationQueries] by element wise product and using array that to take the dot product with [query]
     */
    override fun massPearson(
        timeSeries: F64FlatArray,
        query: F64FlatArray,
        aggregationQueries: F64Array?, // can be supplied to provide aggregation to timeSeries
        reducer: AggregationWithReducerWithArgWithCompare,

        // optional vertical speedup data
        queryDotAggregationProducts: Double?, // can be supplied to prevent double calculations
        queryTimeSeriesSlidingDotProducts: F64FlatArray?,
        queryMean: Double?,
        queryStd: Double?,
        timeSeriesSlidingMeans: F64FlatArray?,
        timeSeriesSlidingStds: F64FlatArray?,
        precalculatedTimeSeriesSlidingMeans: F64FlatArray?, // for calculating sliding mean/std with aggregation faster

        // optional horizontal speedup data
        considerCacheQueryTimeSeriesSlidingDotProductsEmpty: Delegate<Boolean>?,
        cacheQueryTimeSeriesSlidingDotProducts: F64FlatArray?,
        previousQuery: F64FlatArray?,


        // get optional result line
        resultLineCallback: ((F64FlatArray) -> Unit)?,

        // get optional result index
        resultIndexCallback: ((Int) -> Unit)?,
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
            queryTimeSeriesSlidingDotProducts != null -> queryTimeSeriesSlidingDotProducts
            QTCache != null && considerCacheQueryTimeSeriesSlidingDotProductsEmpty != null -> QTCache.also {
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
            else -> slidingDotProducts(Q = Q, T = T, returnArraySlot = Slot.B3)
        }

        // apply aggregation for T on QT
        if (Qa != null && QDotsQa != null) {
            QT += QDotsQa
            QT /= q
        }

        // compute mean and stds
        val mu_Q = queryMean ?: Q.mean()
        val sigma_Q = queryStd
            ?: Q.let {
                it -= mu_Q
                sqrt(it.dot() / it.length)
            }

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

        return P[index]
    }

    /**
     * TODO Lag bounds will mess up qtcache?
     * TODO continuity, so if T has NaNs
     */
    @Suppress("UNUSED_VALUE")
    override fun updateQTCache(
        QTCache: F64FlatArray,
        considerQTCacheEmpty: Delegate<Boolean>,
        Q: F64FlatArray,
        T: F64FlatArray,
        n: Int,
        m: Int,
        previousQuery: F64FlatArray?,
        returnArraySlot: Slot,
    ) {
        require(T.all { it.isNotNaN() }) { "updateQTCache does not (yet) work if it contains NaNs." }
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

    override fun computeSlidingMeanStdWithAgg(
        T: F64FlatArray,
        windowSize: Int,
        aggregations: F64Array,
        precalculatedTMeans: F64FlatArray?,
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
     * Compatible with continuity removed (NaNs in [T]).
     */
    override fun computeSlidingMeanStdWithAgg(
        T: F64FlatArray,
        windowSize: Int,
        aggregations: F64Array,
        destMeanArray: F64FlatArray,
        destStdArray: F64FlatArray,
        precalculatedTMeans: F64FlatArray?,
    ) {
        if (T.any { it.isNaN() }) println("WARNING: slidingMeanAndStdWithAggregation is called with a continuity-removed time series. Results might not be as fast.")

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

            // Tw.dot() == squareBuffer
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

}