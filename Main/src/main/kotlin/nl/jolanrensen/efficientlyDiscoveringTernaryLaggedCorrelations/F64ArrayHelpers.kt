package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations

import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor.asF64Array
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt
import kotlin.random.Random

fun F64FlatArray.zNormalize() {
    val mean = mean()
    val std = standardDeviation()

    this -= mean
    this /= std
}

fun F64FlatArray.zNormalized(): F64FlatArray = copy().apply { zNormalize() }

fun F64FlatArray.standardDeviation(): Double = sqrt(sumOf { (it - mean()).pow(2.0) } / length)

infix fun F64Array.concat(other: F64Array): F64Array = F64Array.concatenate(this, other)
infix fun F64FlatArray.concat(other: F64FlatArray): F64FlatArray = F64Array.concatenate(this, other) as F64FlatArray

inline fun F64FlatArray.sumOf(modifier: (Double) -> Double): Double {
    var sum = 0.0
    for (i in 0 until length) {
        sum += modifier(this[i])
    }
    return sum
}

inline operator fun F64FlatArray.invoke(fill: (Int) -> Double): F64FlatArray {
    for (i in indices)
        this[i] = fill(i)

    return this
}

inline val F64FlatArray.indices: IntRange
    get() = 0 until length

@Suppress("NOTHING_TO_INLINE")
inline fun F64Array.asFlatArray(): F64FlatArray {
    require(nDim == 1)
    return this as F64FlatArray
}

/**
 * Creates a sliced view of this array in O(1) time.
 *
 * @param from the first index of the slice (inclusive).
 * @param to the last index of the slice (exclusive). `-1` is treated as "until the end", otherwise [to] must be
 * strictly greater than [from] (empty arrays are not allowed).
 * @param step indexing step.
 */
@Suppress("NOTHING_TO_INLINE")
inline fun F64FlatArray.sliceFlat(from: Int = 0, to: Int = -1, step: Int = 1): F64FlatArray =
    slice(from = from, to = to, step = step, axis = 0).asFlatArray()

@Suppress("NOTHING_TO_INLINE")
inline fun F64FlatArray.sliceFlat(range: IntRange): F64FlatArray =
    slice(from = range.first, to = range.last + 1, step = 1, axis = 0).asFlatArray()

inline fun F64FlatArray.any(predicate: (Double) -> Boolean): Boolean {
    if (length == 0) return false
    for (i in 0 until length) if (predicate(this[i])) return true
    return false
}

inline fun F64FlatArray.all(predicate: (Double) -> Boolean): Boolean {
    if (length == 0) return true
    for (i in 0 until length) if (!predicate(this[i])) return false
    return true
}

@Suppress("NOTHING_TO_INLINE")
inline fun checkCountOverflow(count: Int): Int {
    if (count < 0) {
        throw ArithmeticException("Count overflow has happened.")
    }
    return count
}

inline fun F64FlatArray.count(predicate: (Double) -> Boolean): Int {
    if (length == 0) return 0
    var count = 0
    for (i in 0 until length) if (predicate(this[i])) checkCountOverflow(++count)
    return count
}

operator fun F64FlatArray.iterator(): Iterator<Double> = object : Iterator<Double> {
    var i = 0
    override fun hasNext(): Boolean = i < length - 1
    override fun next(): Double {
        val res = this@iterator[i]
        i++
        return res
    }
}

fun F64FlatArray.withIndex(): Iterator<IndexedValue<Double>> = iterator().withIndex()

fun F64FlatArray.asIterable(): F64FlatArrayIterable = F64FlatArrayIterable(this)

@JvmInline
value class F64FlatArrayIterable(val array: F64FlatArray) : Iterable<Double> {
    override fun iterator(): Iterator<Double> = array.iterator()
}

fun F64Array.asElementwiseIterable(): Iterable<Double> = flatten().asIterable()

/**
 * Returns the dot product of receiver array with itself
 */
fun F64FlatArray.dot(): Double = this dot this

/**
 * Returns the product of the elements.
 * Optimized for dense arrays.
 */
fun F64FlatArray.product(): Double = if (length == 1) this[0] else reduce(Double::times)

//inline fun F64FlatArray.reduce(operation: (acc: Double, Double) -> Double): Double {
//    if (size == 0) throw UnsupportedOperationException("Empty array can't be reduced.")
//
//    var accumulator = this[0]
//    var i = 1
//    while (i < size) {
//        accumulator = operation(accumulator, this[i++])
//    }
//    return accumulator
//}

/**
 * Performs element^[x] in place for each value of receiver array.
 * @return receiver array itself
 */
fun <T : F64Array> T.powInPlace(x: Double): T = apply { transformInPlace { it.pow(x) } }

/**
 * Performs element^[x] in copy of receiver array for each value.
 * @return new array
 */
@Suppress("UNCHECKED_CAST")
fun <T : F64Array> T.pow(x: Double): T = copy().apply { powInPlace(x) } as T


/**
 * Performs sqrt(element) in place for each value of receiver array.
 * @return receiver array itself
 */
fun <T : F64Array> T.sqrtInPlace(): T = apply { transformInPlace { sqrt(it) } }

/**
 * Performs sqrt(element) in copy of receiver array for each value.
 * @return new array
 */
@Suppress("UNCHECKED_CAST")
fun <T : F64Array> T.sqrt(): T = copy().apply { sqrtInPlace() } as T

/**
 * Creates an [F64FlatArray] of given [size] filled with random doubles.
 */
fun randomF64FlatArray(size: Int, from: Double? = null, until: Double? = null): F64FlatArray =
    F64FlatArray(size) {
        when {
            from == null && until == null -> Random.nextDouble()
            from == null && until != null -> Random.nextDouble(until = until)
            from != null && until == null -> throw IllegalArgumentException()
            else -> Random.nextDouble(from = from!!, until = until!!)
        }
    }

/**
 * Creates a new [F64FlatArray] filled with 0 of size [size].
 */
operator fun F64FlatArray.Companion.invoke(size: Int): F64FlatArray = DoubleArray(size = size).asF64Array()

/**
 * Creates a new [F64FlatArray] of size [size] which can be filled up using [block].
 */
operator fun F64FlatArray.Companion.invoke(size: Int, block: (Int) -> Double): F64FlatArray =
    invoke(size).apply {
        for (i in 0 until size) {
            this[i] = block(i)
        }
    }

/**
 * Creates a new [F64FlatArray] filled with [init] of size [size].
 */
fun F64FlatArray.Companion.full(size: Int, init: Double): F64FlatArray = invoke(size).apply { fill(init) }

/**
 * Provides a view of receiver array and tries to cast result to an [F64FlatArray].
 * @throws IllegalArgumentException if this fails.
 */
@Suppress("NOTHING_TO_INLINE")
inline fun F64Array.viewAsFlat(index: Int, axis: Int = 0): F64FlatArray = view(index = index, axis = axis).asFlatArray()

/**
 * Reverses receiver array in place.
 */
fun F64FlatArray.reverse() {
    var i = length - 1
    var j = 0
    while (i > j) {
        val temp = this[i]
        this[i] = this[j]
        this[j] = temp
        i--
        j++
    }
}

/**
 * Returns a reversed copy of receiver array.
 */
fun F64FlatArray.reversed(): F64FlatArray = copy().also { it.reverse() }

/**
 * Creates a vector from given elements.
 * same as [F64Array.of]
 */
@Deprecated("", ReplaceWith("f64FlatArrayOf"))
fun F64FlatArray.Companion.of(first: Double, vararg rest: Double): F64FlatArray {
    val data = DoubleArray(rest.size + 1)
    data[0] = first
    System.arraycopy(rest, 0, data, 1, rest.size)
    return data.asF64Array()
}

/**
 * Creates a vector from given elements.
 * same as [F64Array.of]
 */
fun f64FlatArrayOf(first: Double, vararg rest: Double): F64FlatArray {
    val data = DoubleArray(rest.size + 1)
    data[0] = first
    System.arraycopy(rest, 0, data, 1, rest.size)
    return data.asF64Array()
}

/**
 * Creates a vector from given elements.
 * same as [F64Array.of]
 */
fun f64FlatArrayOf(first: Number, vararg rest: Number): F64FlatArray {
    val data = DoubleArray(rest.size + 1) {
        if (it == 0) first.toDouble()
        else rest[it - 1].toDouble()
    }
    return data.asF64Array()
}

inline fun <T> F64FlatArray.mapWindowed(
    windowSize: Int,
    block: (query: F64FlatArray) -> T,
): List<T> {
    require(windowSize <= length) { "Window size must be <= size." }
    require(windowSize > 0) {
        if (windowSize != 1) "size $windowSize must be greater than zero."
        else "size $windowSize must be greater than zero."
    }

    return List(length - windowSize + 1) {
        val window = sliceFlat(from = it, to = it + windowSize)
        block(window)
    }
}

inline fun <T> F64FlatArray.mapWindowedWithIndex(
    windowSize: Int,
    block: (i: Int, query: F64FlatArray) -> T,
): List<T> {
    var i = 0
    return mapWindowed(windowSize) {
        block(i++, it)
    }
}

/**
 * Executes [block] for all placements of a window of size [windowSize].
 */
inline fun F64FlatArray.forEachWindow(
    windowSize: Int,
    block: (query: F64FlatArray) -> Unit,
) {
    require(windowSize <= length) { "Window size must be <= size." }
    require(windowSize > 0) {
        if (windowSize != 1) "size $windowSize must be greater than zero."
        else "size $windowSize must be greater than zero."
    }

    for (i in 0..length - windowSize) {
        val window = sliceFlat(from = i, to = i + windowSize)
        block(window)
    }
}

/**
 * Executes [block] for all placements of a window of size [windowSize].
 * Includes index [i].
 */
inline fun F64FlatArray.forEachWindowWithIndex(
    windowSize: Int,
    block: (i: Int, query: F64FlatArray) -> Unit,
) {
    var i = 0
    forEachWindow(windowSize) {
        block(i++, it)
    }
}

fun F64FlatArray.toList(): List<Double> = toDoubleArray().toList()

inline fun F64FlatArray.transformInPlaceWithIndex(crossinline op: (Int, Double) -> Double) {
    var i = 0
    transformInPlace { op(i++, it) }
}

inline fun F64FlatArray.transformWithIndex(crossinline op: (Int, Double) -> Double): F64FlatArray =
    copy().apply { transformInPlaceWithIndex(op) }

fun F64Array.transposed(): F64Array = copy().apply { transposeInPlace() }

fun F64Array.transposeInPlace() {
    require(shape.size == 1 || shape.size == 2)

    when (shape.size) {
        1 -> asFlatArray().reverse()
        2 -> {
            for (row in 0 until length) {
                for (col in row until V[row].length) {
                    val data = this[row, col]
                    this[row, col] = this[col, row]
                    this[col, row] = data
                }
            }
        }
    }
}

fun getSlopes(array: F64FlatArray) = F64FlatArray(array.length - 1) { array[it + 1] - array[it] }

/**
 * Returns a new array with slopes and curve points denoted by NaN
 */
fun F64FlatArray.getSlopesWithCurvePoints(): F64FlatArray {
    val list = mutableListOf<Double>()
    val slopes = getSlopes(this)

    var prev: Double? = null
    for (item in slopes) {
        if (prev != null && abs(item - prev) > 0.001) list += Double.NaN
        list += item
        prev = item
    }

    return list.toDoubleArray().asF64Array()
}

fun F64FlatArray.encodeToString(): String = toDoubleArray().asIterable().joinToString(";")

fun F64FlatArray.Companion.decodeFromString(string: String): F64FlatArray =
    string.split(';')
        .map { it.toDouble() }
        .toDoubleArray()
        .asF64Array()

fun F64FlatArray.equalsRoughly(other: F64FlatArray, precision: Double = 0.0001): Boolean = when {
    this === other -> true
    length != other.length -> false
    else -> (0 until length).all {
        this[it].closeTo(other[it], precision)
    }
}