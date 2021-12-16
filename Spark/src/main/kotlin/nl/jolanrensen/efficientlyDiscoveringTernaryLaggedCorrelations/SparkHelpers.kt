@file:Suppress("RemoveExplicitTypeArguments")

package nl.jolanrensen.thesis

import nl.jolanrensen.scalaTuplesInKotlin.first
import nl.jolanrensen.scalaTuplesInKotlin.second
import nl.jolanrensen.scalaTuplesInKotlin.u
import org.apache.spark.api.java.JavaPairRDD
import org.apache.spark.api.java.JavaRDDLike
import org.apache.spark.api.java.function.*
import org.apache.spark.broadcast.Broadcast
import org.apache.spark.sql.Column
import org.apache.spark.sql.Dataset
import org.jetbrains.kotlinx.spark.api.Arity2
import org.jetbrains.kotlinx.spark.api.encoder
import org.jetbrains.kotlinx.spark.api.forEach
import org.jetbrains.kotlinx.spark.api.map
import scala.Tuple2
import kotlin.random.Random
import kotlin.reflect.KProperty
import kotlin.reflect.KProperty1

@Deprecated("This breaks if forEach is done parallel")
inline fun <reified T> Dataset<T>.forEachIndexed(noinline action: (index: Int, T) -> Unit) {
    var index = 0
    forEach {
        action(index++, it)
    }
}

/** Performs for each in partitions on receiver and for each on iterator. */
inline fun <T> Dataset<T>.forEachInPartitions(noinline func: (T) -> Unit): Unit =
    foreachPartition(ForeachPartitionFunction { it.forEach(func) })

inline fun <T> Dataset<T>.forEachPartition(noinline func: (Iterator<T>) -> Unit): Unit =
    foreachPartition(func)

/** Maps receiver in partitions and maps iterator. */
inline fun <T, reified R> Dataset<T>.mapInPartitions(noinline func: (T) -> R): Dataset<R> =
    mapPartitions(MapPartitionsFunction { it.map(func) }, encoder<R>())


operator fun <T> Broadcast<T>.getValue(thisRef: Any?, property: KProperty<*>): T = value()


//inline fun <reified T1, T2> Dataset<Tuple2<T1, T2>>.takeKeys(): Dataset<T1> = map { it.first }
@JvmName("takeKeysInPartitionsTuple2")
inline fun <reified T1, T2> Dataset<Tuple2<T1, T2>>.takeKeysInPartitions(): Dataset<T1> =
    mapInPartitions { it.first }

inline fun <reified T1, T2> Dataset<Pair<T1, T2>>.takeKeysInPartitions(): Dataset<T1> =
    mapInPartitions { it.first }

@JvmName("takeValuesInPartitionsTuple2")
inline fun <T1, reified T2> Dataset<Tuple2<T1, T2>>.takeValuesInPartitions(): Dataset<T2> =
    mapInPartitions { it.second }

inline fun <T1, reified T2> Dataset<Pair<T1, T2>>.takeValuesInPartitions(): Dataset<T2> =
    mapInPartitions { it.second }


// TODO just experimenting
class KDataset<T>(dataset: Dataset<T>) : Dataset<T>(dataset.queryExecution(), dataset.encoder()) {
    fun reduce(func: (T, T) -> T): T = reduce(ReduceFunction(func))
}

fun <T> Dataset<T>.toKotlin(): KDataset<T> = KDataset(this)

fun <T> emptyIterator(): Iterator<T> = emptyList<Nothing>().iterator()

fun <T> iteratorOf(arg: T): Iterator<T> = object : Iterator<T> {
    var taken = false
    override fun hasNext(): Boolean = !taken

    override fun next(): T {
        if (taken) {
            throw IndexOutOfBoundsException()
        } else {
            taken = true
            return arg
        }
    }
}

fun <T> iteratorOf(vararg args: T): Iterator<T> = if (args.isEmpty()) emptyIterator<T>() else args.iterator()

/** Allows for `YourClass::property.col()`. */
fun KProperty1<*, *>.col(): Column = org.apache.spark.sql.functions.col(name)

@Suppress("RedundantSamConstructor")
inline fun <reified T, This : JavaRDDLike<T, This>, reified K2, reified V2> JavaRDDLike<T, This>.mapPartitionsToPairK(
    noinline func: (t: Iterator<T>) -> Iterator<Tuple2<K2, V2>>,
): JavaPairRDD<K2, V2> = mapPartitionsToPair(PairFlatMapFunction(func)) as JavaPairRDD<K2, V2>

@Suppress("RedundantSamConstructor")
inline fun <reified K, reified V, reified U> JavaPairRDD<K, V>.aggregateByKey(
    zeroValue: U,
    noinline seqFunction: (U, V) -> U,
    noinline combineFunction: (U, U) -> U,
): JavaPairRDD<K, U> = aggregateByKey<U>(
    zeroValue,
    Function2(seqFunction),
    Function2(combineFunction),
) as JavaPairRDD<K, U>

class DatasetSequence<T>(dataset: Dataset<T>) : Dataset<T>(dataset.queryExecution(), dataset.encoder()), Sequence<T> {
    override fun iterator(): Iterator<T> = toLocalIterator()
}

class DatasetIterable<T>(dataset: Dataset<T>) : Dataset<T>(dataset.queryExecution(), dataset.encoder()), Iterable<T> {
    override fun iterator(): Iterator<T> = toLocalIterator()
}

inline fun <reified T> Dataset<T>.asSequence(): DatasetSequence<T> = DatasetSequence(this)
inline fun <reified T> Dataset<T>.asIterable(): DatasetIterable<T> = DatasetIterable(this)

inline fun <reified T> Iterator<T>.asIterable(): Iterable<T> = Iterable { this }

inline fun <reified T> Dataset<T>.randomized(): Dataset<T> =
    mapInPartitions { it u Random.nextInt() }
        .orderBy(org.jetbrains.kotlinx.spark.api.col(Tuple2<*, *>::_2))
        .mapInPartitions { it.first }

@JvmName("sortByKeyTuple2")
fun <T1, T2> Dataset<Tuple2<T1, T2>>.sortByKey(): Dataset<Tuple2<T1, T2>> = sort("_1")

@JvmName("sortByValueTuple2")
fun <T1, T2> Dataset<Tuple2<T1, T2>>.sortByValue(): Dataset<Tuple2<T1, T2>> = sort("_2")

@JvmName("sortByKeyArity2")
fun <T1, T2> Dataset<Arity2<T1, T2>>.sortByKey(): Dataset<Arity2<T1, T2>> = sort("_1")

@JvmName("sortByValueArity2")
fun <T1, T2> Dataset<Arity2<T1, T2>>.sortByValue(): Dataset<Arity2<T1, T2>> = sort("_2")

fun <T1, T2> Dataset<Pair<T1, T2>>.sortByKey(): Dataset<Pair<T1, T2>> = sort("first")

fun <T1, T2> Dataset<Pair<T1, T2>>.sortByValue(): Dataset<Pair<T1, T2>> = sort("second")


