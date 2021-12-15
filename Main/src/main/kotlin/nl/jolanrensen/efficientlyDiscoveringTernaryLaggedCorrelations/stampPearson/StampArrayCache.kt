package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.invoke
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.nextPowerOf2
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.sliceFlat
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampArrayCache.Slot.*
import org.jetbrains.bio.viktor.F64FlatArray
import java.io.Serializable

/**
 *
 * Array cache used by [StampPearson], [StampPearson3ts], [StampPearson3tsWithSkipping]
 * This makes sure that no new memory needs to be allocated for a second stamp run.
 * It does make it so that one stamp instance cannot run in parallel.
 *
 * @param maxSize usually largest time series size
 *
 */
class StampArrayCache(
    private val maxSize: Int = 100,
    private val slotManager: SlotManager = SlotsWithReuse,
) : Serializable {

    enum class Slot : Serializable {
        A, // used in slidingDotProducts
        B, // used in slidingDotProducts
        B1, // queryCTimeSeriesASlidingDotProducts
        B2, // queryCTimeSeriesBSlidingDotProducts
        B3, // QT
        B4, // TQaSumsSlidingDots
        C, // cacheQueryCTimeSeriesASlidingDotProducts
        D, // cacheQueryCTimeSeriesBSlidingDotProducts
        E, // cacheQueryBTimeSeriesASlidingDotProducts
        F, // timeSeriesASlidingMeans
        G, // timeSeriesASlidingStds
        H, // timeSeriesASlidingMeansAggQC
        I, // timeSeriesASlidingStdsAggQC
        J, // timeSeriesBSlidingMeans
        K, // timeSeriesBSlidingMeans
        L, // cacheQueryBTimeSeriesASlidingDotProducts
        M, // queryCStd
        N, // agg(QB, QC)
        O, // mus_T
        P, // sigmas_T
        Q, // QaSums
        R, // QDotsQa
        S, // timeSeriesCSlidingMeans
        T, // timeSeriesCSlidingStds
        U, // query = QC
    }

    sealed interface SlotManager : Serializable {
        operator fun invoke(slot: Slot): Int
    }

    /**
     * Has some slots overlapping to save space.
     * Default SlotManager.
     */
    object SlotsWithReuse : SlotManager {
        override fun invoke(slot: Slot): Int = when (slot) {
            A -> 0
            B -> 1
            B1 -> 17
            B2 -> 18
            B3 -> 19
            B4 -> 20
            C -> 2
            D -> 3
            E -> 4
            F -> 5
            G -> 6
            H -> 7
            I -> 8
            J -> 9
            K -> 10
            L -> 11
            M -> 0
            N -> 12
            O -> 0
            P -> 1
            Q -> 13
            R -> 0
            S -> 14
            T -> 15
            U -> 16
        }
    }

    /**
     * All slots are unique, better for debugging at cost of space.
     */
    object SlotsWithoutReuse : SlotManager {
        override fun invoke(slot: Slot): Int = when (slot) {
            A -> 0
            B -> 1
            B1 -> 21
            B2 -> 22
            B3 -> 23
            B4 -> 24
            C -> 2
            D -> 3
            E -> 4
            F -> 5
            G -> 6
            H -> 7
            I -> 8
            J -> 9
            K -> 10
            L -> 11
            M -> 12
            N -> 13
            O -> 14
            P -> 15
            Q -> 16
            R -> 17
            S -> 18
            T -> 19
            U -> 20
        }
    }

    /**
     * Array cache is not used.
     */
    object NoSlots : SlotManager {
        override fun invoke(slot: Slot): Int = -1
    }

    private val arrayCache: MutableMap<Int, F64FlatArray> =
        when (slotManager) {
            SlotsWithReuse, SlotsWithoutReuse ->
                Slot.values()
                    .map {
                        Pair(slotManager(it), it in listOf(A, B, B1, B2, B3, B4))
                    }
                    .associateWith { (_, isSlidingDots) ->
                        F64FlatArray(
                            if (isSlidingDots) 2 * nextPowerOf2(maxSize)
                            else maxSize
                        )
                    }.mapKeys { (key, _) -> key.first }
                    .toMutableMap()

            NoSlots -> mutableMapOf()
        }


    operator fun get(size: Int, slot: Slot): F64FlatArray = when {

        slotManager == NoSlots ->
            F64FlatArray(size)

        slotManager(slot) !in arrayCache ->
            F64FlatArray(size)
                .also { arrayCache[slotManager(slot)] = it }

        arrayCache[slotManager(slot)]!!.length == size ->
            arrayCache[slotManager(slot)]!!

        arrayCache[slotManager(slot)]!!.length > size ->
            arrayCache[slotManager(slot)]!!
                .sliceFlat(0, size)

        else ->
            F64FlatArray(size).also {
                arrayCache[slotManager(slot)] = it
                println("had to enlarge slot $slot to size $size")
            }
    }

    override fun toString(): String = arrayCache.toString()
}