package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations

import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor.minus
import org.jetbrains.bio.viktor.times
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt


fun nextPowerOf2(a: Int): Int {
    var b = 1
    while (b < a) {
        b = b shl 1
    }
    return b
}

inline fun <T, R> Iterable<T>.windowedWithIndex(
    size: Int,
    step: Int = 1,
    partialWindows: Boolean = false,
    crossinline transform: (Int, List<T>) -> R,
): List<R> {
    var i = 0
    return windowed(size = size, step = step, partialWindows = partialWindows) { window: List<T> ->
        transform(i++, window)
    }
}

fun euclideanToPearsonCorrelations(series: F64FlatArray, windowSize: Int): F64FlatArray =
    (1.0 - (1.0 / (2.0 * windowSize)) * (series * series)).asFlatArray()


inline fun Double.isNotNaN(): Boolean = !this.isNaN()
inline fun Float.isNotNaN(): Boolean = !this.isNaN()

inline fun Double.coerceNotNaN(replaceWith: Double = 0.0): Double = if (isNaN()) replaceWith else this
inline fun Float.coerceNotNaN(replaceWith: Float = 0f): Float = if (isNaN()) replaceWith else this

inline fun Double.coerceNotInf(replaceWith: Double = 0.0): Double = if (isInfinite()) replaceWith else this
inline fun Float.coerceNotInf(replaceWith: Float = 0f): Float = if (isInfinite()) replaceWith else this

fun Double.toString(noDecimals: Int? = null, totalLength: Int? = null): String =
    let {
        when {
            noDecimals == null -> toString()
            noDecimals > 0 -> "%.${noDecimals}f".format(this)
            else -> throw IllegalArgumentException("noDecimals must be > 0.")
        }
    }.let {
        when {
            totalLength == null -> it
            else -> it.padStart(length = totalLength)
        }
    }

/**
 * Splits the original collection into pair of lists,
 * where *first* list contains elements for which [predicate] yielded `true`,
 * while *second* list contains elements for which [predicate] yielded `false`.
 *
 * @sample samples.collections.Iterables.Operations.partition
 */
inline fun <T> Iterable<T>.partitionIndexed(predicate: (Int, T) -> Boolean): Pair<List<T>, List<T>> {
    val first = ArrayList<T>()
    val second = ArrayList<T>()
    for ((i, element) in this.withIndex()) {
        if (predicate(i, element)) {
            first.add(element)
        } else {
            second.add(element)
        }
    }
    return Pair(first, second)
}

inline fun Double.closeTo(other: Double, precision: Double): Boolean = abs(this - other) <= precision
inline fun Double.notCloseTo(other: Double, precision: Double): Boolean = abs(this - other) > precision

inline fun <T> Iterable<T>.collectionSizeOrDefault(default: Int): Int =
    if (this is Collection<*>) this.size else default

public fun <T, R, S> Iterable<Triple<T, R, S>>.unzip(): Triple<List<T>, List<R>, List<S>> {
    val expectedSize = collectionSizeOrDefault(10)
    val listT = ArrayList<T>(expectedSize)
    val listR = ArrayList<R>(expectedSize)
    val listS = ArrayList<S>(expectedSize)
    for (triple in this) {
        listT.add(triple.first)
        listR.add(triple.second)
        listS.add(triple.third)
    }
    return Triple(listT, listR, listS)
}

inline infix fun <T : Comparable<T>> ClosedRange<T>.overlaps(other: ClosedRange<T>): Boolean =
    this.endInclusive >= other.start && this.start <= other.endInclusive

class PolyRegression(
    private val a: Double,
    private val b: Double,
    private val c: Double,
    x: F64FlatArray,
    y: F64FlatArray,
) {

    init {
        require(x.length == y.length)
    }

    val n: Int = x.length

    val r2: Double = (
            (n * (x * y).sum() - x.sum() * y.sum()
                    ) / (
                    sqrt(n * x.pow(2.0).sum() - x.sum().pow(2.0)) *
                            sqrt(n * y.pow(2.0).sum() - y.sum().pow(2.0))
                    )
            ).pow(2.0)


    fun function(x: Double) = a + b * x + c * x * x

    operator fun invoke(x: Double) = function(x)
}

fun polyRegression(x: F64FlatArray, y: F64FlatArray): PolyRegression {
    val xm = x.mean()
    val ym = y.mean()

//    val x2m = x.map { it * it }.average()
    val x2m = x.pow(2.0).mean()

//    val x3m = x.map { it * it * it }.average()
    val x3m = x.pow(3.0).mean()

//    val x4m = x.map { it * it * it * it }.average()
    val x4m = x.pow(4.0).mean()

//    val xym = x.zip(y).map { it.first * it.second }.average()
    val xym = (x * y).mean()

//    val x2ym = x.zip(y).map { it.first * it.first * it.second }.average()
    val x2ym = (x.pow(2.0) * y).mean()

    val sxx = x2m - xm * xm
    val sxy = xym - xm * ym
    val sxx2 = x3m - xm * x2m
    val sx2x2 = x4m - x2m * x2m
    val sx2y = x2ym - x2m * ym

    val b = (sxy * sx2x2 - sx2y * sxx2) / (sxx * sx2x2 - sxx2 * sxx2)
    val c = (sx2y * sxx - sxy * sxx2) / (sxx * sx2x2 - sxx2 * sxx2)
    val a = ym - b * xm - c * x2m

    val result = PolyRegression(a, b, c, x, y)

    println("y = $a + ${b}x + ${c}x^2\n")
    println(" Input  Approximation")
    println(" x   y     y1")
    for (i in x.indices) {
        val xi = x[i]
        val yi = y[i]

        println("$xi $yi  ${result(xi)}")
    }

    return result
}