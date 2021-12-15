package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.asElementwiseIterable
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import java.io.Serializable


sealed interface Aggregation : Serializable {

    // To be implemented

    /** Aggregates given [F64Array] to a single [Double] value. */
    fun aggregate(array: F64Array): Double

    /** Aggregates given [DoubleArray] to a single [Double] value. */
    fun aggregate(array: DoubleArray): Double

    // Helper

    /** Aggregates given [F64Array] to a single [Double] value. */
    operator fun invoke(array: F64Array): Double = aggregate(array)

    /** Aggregates given [DoubleArray] to a single [Double] value. */
    operator fun invoke(array: DoubleArray): Double = aggregate(array)
}


sealed interface AggregationWithReducer : Aggregation {

    // To be implemented

    /** Reduce [Double]s [a] and [b] in the same way the [aggregate] method does. */
    fun reduce(a: Double, b: Double): Double

    /** The "starting" value of the reducer. For example: infinity for [MIN], 1.0 for [PRODUCT], 0.0 for [SUM]. */
    val initializer: Double


    // Helper functions

    /** Replaces all [Double.NaN] values in [array] with [initializer] value. */
    fun replaceNaN(array: F64Array): Unit = array.transformInPlace(::coerceNotNaN)

    /** Replaces all [Double.NaN] values in [array] with [initializer] value. */
    fun replaceNaN(array: DoubleArray) {
        for (i in array.indices) {
            array[i] = coerceNotNaN(array[i])
        }
    }

    fun coerceNotNaN(value: Double): Double = when {
        value.isNaN() -> initializer
        else -> value
    }

    /** Reduce [Double]s [a] and [b] in the same way the [aggregate] method does. */
    operator fun invoke(a: Double, b: Double): Double = reduce(a, b)

    /** Reduce [Double]s [a], [b], and [c] in the same way the [aggregate] method does. */
    fun reduce(a: Double, b: Double, c: Double): Double = reduce(a, reduce(b, c))

    /** Reduce [Double]s [a], [b], and [c] in the same way the [aggregate] method does. */
    operator fun invoke(a: Double, b: Double, c: Double): Double = reduce(a, reduce(b, c))

    /** Reduce [Double]s [a], [b], [c], and [d] in the same way the [aggregate] method does. */
    fun reduce(a: Double, b: Double, c: Double, d: Double) = reduce(reduce(a, b), reduce(c, d))

    /** Reduce [Double]s [a], [b], [c], and [d] in the same way the [aggregate] method does. */
    operator fun invoke(a: Double, b: Double, c: Double, d: Double) = reduce(reduce(a, b), reduce(c, d))
}

sealed interface AggregationWithArg : Aggregation {

    // To be implemented

    /** Returns the index of where the result of [aggregate] is found in [array]. */
    fun argFunction(array: F64Array): Int
}

sealed interface AggregationWithCompare {

    // To be implemented
    val comparator: Comparator<Double>

    fun firstIsBetterThanSecond(first: Double, second: Double): Boolean = comparator.compare(first, second) > 0
}

sealed interface AggregationWithReducerWithArgWithCompare :
    AggregationWithReducer,
    AggregationWithArg,
    AggregationWithCompare


// Implementations

/** Can get the mean/average of an [F64Array]. */
object MEAN : Aggregation {
    override fun aggregate(array: F64Array): Double = array.mean()
    override fun aggregate(array: DoubleArray): Double = array.average()
}

/** Can get the unbiased standard deviation of an [F64Array]. */
object UNBIASED_SD : Aggregation {
    override fun aggregate(array: F64Array): Double = array.sd()
    override fun aggregate(array: DoubleArray): Double = array.asF64Array().sd()
}

/** Has a product reducer and can get the overall product of an [F64Array]. */
object PRODUCT : AggregationWithReducer {
    override fun reduce(a: Double, b: Double): Double = a * b
    override val initializer: Double = 1.0
    override fun aggregate(array: F64Array): Double = array.asElementwiseIterable().reduce(::reduce)
    override fun aggregate(array: DoubleArray): Double = array.reduce(PRODUCT::reduce)
}

/** Has a sum reducer and can get the overall sum of an [F64Array]. */
object SUM : AggregationWithReducer {
    override fun aggregate(array: F64Array): Double = array.sum()
    override fun aggregate(array: DoubleArray): Double = array.sum()
    override fun reduce(a: Double, b: Double): Double = a + b
    override val initializer: Double = 0.0
}

/** Has a [minOf] reducer, can get the overall min and the index of that value in an [F64Array]. */
object MIN : AggregationWithReducerWithArgWithCompare {
    override fun aggregate(array: F64Array): Double = array.min()
    override fun aggregate(array: DoubleArray): Double = array.minOrNull() ?: Double.POSITIVE_INFINITY

    override fun argFunction(array: F64Array): Int = array.argMin()

    override val comparator: Comparator<Double> = Comparator { a, b -> -a.compareTo(b) }

    override fun reduce(a: Double, b: Double): Double = minOf(a, b)
    override val initializer: Double = Double.POSITIVE_INFINITY
}

/** Has a [maxOf] reducer, can get the overall max and the index of that value in an [F64Array]. */
object MAX : AggregationWithReducerWithArgWithCompare {
    override fun aggregate(array: F64Array): Double = array.max()
    override fun aggregate(array: DoubleArray): Double = array.maxOrNull() ?: Double.NEGATIVE_INFINITY

    override fun argFunction(array: F64Array): Int = array.argMax()

    override val comparator: Comparator<Double> = Comparator { a, b -> a.compareTo(b) }

    override fun reduce(a: Double, b: Double): Double = maxOf(a, b)
    override val initializer: Double = Double.NEGATIVE_INFINITY
}
