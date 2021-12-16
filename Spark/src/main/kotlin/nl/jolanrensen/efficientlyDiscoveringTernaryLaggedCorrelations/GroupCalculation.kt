package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations

import nl.jolanrensen.scalaTuplesInKotlin.tupleOf
import nl.jolanrensen.thesis.aggregateByKey
import nl.jolanrensen.thesis.mapPartitionsToPairK
import org.apache.spark.api.java.JavaRDD
import org.apache.spark.api.java.JavaSparkContext
import org.apache.spark.api.java.function.Function2
import org.apache.spark.sql.Dataset
import org.jetbrains.kotlinx.spark.api.*
import scala.Tuple2
import java.util.*
import java.util.stream.Collectors
import java.util.stream.IntStream
import kotlin.math.pow

object GroupCalculation {
    /**
     * Simple method to give each index of x dimensions a unique number.
     *
     * @param indexTuple int[] (or tuple) of size x with all values < listSize. The index for which to return the number
     * @param listSize   the size of the list, aka the max width, height etc of the table
     * @return the unique number for this `indexTuple`
     */
    private fun getTupleValue(indexTuple: IntArray, listSize: Int): Int {
        var result = 0
        for (i in indexTuple.indices) {
            result += indexTuple[i] * listSize.toDouble().pow(i.toDouble()).toInt()
        }
        return result
    }

    /**
     * To make sure that every tuple is only picked once, this method returns true only if the indices are in the right
     * corner of the matrix. This works for any number of dimensions > 1. Here is an example for 2-D:
     *
     *
     * -  0  1  2  3  4  5  6  7  8  9
     * --------------------------------
     * 0| x  ✓  ✓  ✓  ✓  ✓  ✓  ✓  ✓  ✓
     * 1| x  x  ✓  ✓  ✓  ✓  ✓  ✓  ✓  ✓
     * 2| x  x  x  ✓  ✓  ✓  ✓  ✓  ✓  ✓
     * 3| x  x  x  x  ✓  ✓  ✓  ✓  ✓  ✓
     * 4| x  x  x  x  x  ✓  ✓  ✓  ✓  ✓
     * 5| x  x  x  x  x  x  ✓  ✓  ✓  ✓
     * 6| x  x  x  x  x  x  x  ✓  ✓  ✓
     * 7| x  x  x  x  x  x  x  x  ✓  ✓
     * 8| x  x  x  x  x  x  x  x  x  ✓
     * 9| x  x  x  x  x  x  x  x  x  x
     *
     * @param indexTuple a tuple of indices in the form of an int[]
     * @return true if this tuple is in the right corner and should be included
     */
    private fun isValidIndexTuple(indexTuple: IntArray): Boolean {
        // x - y > 0  // 2d
        // (x - y) > 0 && (x - z) > 0 && (y - z) > 0;  // 3d
        // (x - y) > 0 && (x - z) > 0 && (x - a) > 0 && (y - z) > 0  && (y - a) > 0 && (z - a) > 0  // 4d
        require(indexTuple.size >= 2) { "not a tuple" }
        for (i in 0 until indexTuple.size - 1) {
            for (j in i + 1 until indexTuple.size) {
                if (indexTuple[i] - indexTuple[j] <= 0) return false
            }
        }
        return true
    }

    /**
     * Recursive method that for `skipDimension` loops over all the other dimensions and stores
     * in `keyValues` all results from `getTupleValue` as key and `value` as value.
     * In the end, `keyValues` will have, for each key in the table below, a value for the key's column, row etc.
     *
     *
     * This is an example for 2D. The letters will be int indices as well (a = 0, b = 1, ..., `listSize`), but help for clarification.
     * The numbers we don't want are filtered out using `isValidIndexTuple`.
     * The actual value of the number in the table comes from `getTupleValue`.
     *
     *
     *
     *
     * -  a  b  c  d  e  f  g  h  i  j
     * --------------------------------
     * a| -  1  2  3  4  5  6  7  8  9
     * b| -  - 12 13 14 15 16 17 18 19
     * c| -  -  - 23 24 25 26 27 28 29
     * d| -  -  -  - 34 35 36 37 38 39
     * e| -  -  -  -  - 45 46 47 48 49
     * f| -  -  -  -  -  - 56 57 58 59
     * g| -  -  -  -  -  -  - 67 68 69
     * h| -  -  -  -  -  -  -  - 78 79
     * i| -  -  -  -  -  -  -  -  - 89
     * j| -  -  -  -  -  -  -  -  -  -
     *
     *
     * @param groupSize        the size of index tuples to form
     * @param keyValues        the ArrayList in which to store the key-value results
     * @param value            the current index to work from (can be seen as a letter in the table above)
     * @param listSize         the size of the list to make
     * @param currentDimension the indicator for which dimension we're currently calculating for (and how deep in the recursion we are)
     * @param indexTuple       the array (or tuple) in which to store the current indices
     * @param skipDimension    the current dimension that will have a set value `value` while looping over the other dimensions
     */
    private fun addTuples(
        groupSize: Int,
        keyValues: ArrayList<Tuple2<Int, Int>>,
        value: Int,
        listSize: Int,
        currentDimension: Int,
        indexTuple: IntArray,
        skipDimension: Int,
    ) {
        if (currentDimension >= groupSize) {  // base case
            if (isValidIndexTuple(indexTuple)) {
                keyValues += tupleOf(
                    first = getTupleValue(indexTuple, listSize),
                    second = value
                )
            }
            return
        }
        if (currentDimension == skipDimension) {
            indexTuple[currentDimension] = value
            addTuples(
                groupSize = groupSize,
                keyValues = keyValues,
                value = value,
                listSize = listSize,
                currentDimension = currentDimension + 1,
                indexTuple = indexTuple,
                skipDimension = skipDimension
            )
        } else {
            for (i in 0 until listSize) {
                indexTuple[currentDimension] = i
                addTuples(
                    groupSize = groupSize,
                    keyValues = keyValues,
                    value = value,
                    listSize = listSize,
                    currentDimension = currentDimension + 1,
                    indexTuple = indexTuple,
                    skipDimension = skipDimension
                )
            }
        }
    }

    /**
     * Simple helper function to get a filled int[] of specified size and fillValue
     *
     * @param size      the size of the int array
     * @param fillValue the value with which to fill up the array
     * @return the array
     */
    private fun getFilledIntArray(size: Int, fillValue: Int): IntArray {
        val array = IntArray(size)
        Arrays.fill(array, fillValue)
        return array
    }

    /**
     * Get all the possible, unique, non repeating groups (of size `groupSize`) of indices for a list of
     * size `listSize`.
     *
     *
     * The workload is evenly distributed by
     *
     * @param listSize  the size of the list for which to calculate the indices
     * @param groupSize the size of a group of indices
     * @return all the possible, unique non repeating groups of indices
     */
    fun KSparkSession.getAllPossibleGroups(
        listSize: Int,
        groupSize: Int,
    ): Dataset<IntArray> = spark.getAllPossibleGroups(listSize = listSize, groupSize = groupSize)

    @JvmName("getAllPossibleGroupsJava")
    fun getAllPossibleGroups(
        spark: SparkSession,
        listSize: Int,
        groupSize: Int,
    ): Dataset<IntArray> = spark.getAllPossibleGroups(listSize = listSize, groupSize = groupSize)

    /**
     * Get all the possible, unique, non repeating groups (of size `groupSize`) of indices for a list of
     * size `listSize`.
     *
     *
     * The workload is evenly distributed by
     *
     * @param listSize  the size of the list for which to calculate the indices
     * @param groupSize the size of a group of indices
     * @return all the possible, unique non repeating groups of indices
     */
    fun SparkSession.getAllPossibleGroups(
        listSize: Int,
        groupSize: Int,
    ): Dataset<IntArray> {
        val sc = JavaSparkContext(sparkContext)

        val broadcastSize = broadcast(listSize)
        val broadcastGroupSize = broadcast(groupSize)
        val indices = sc.parallelize(
            IntStream.range(0, listSize)
                .boxed()
                .collect(Collectors.toList())
        )

        // for a groupsize of 1, no pairing up is needed, so just return the indices converted to int[]'s
        if (groupSize == 1) {
            return createDataset(
                indices
                    .mapPartitions {
                        val list = ArrayList<IntArray>()
                        while (it.hasNext()) {
                            list += intArrayOf(it.next()!!)
                        }
                        list.iterator()
                    }.rdd(),
                encoder(),
            )
        }
        @Suppress("RedundantSamConstructor")
        val keys = indices
            .mapPartitionsToPairK<Int, JavaRDD<Int>, Int, Int> { // this converts all indices of stocks to (number in table, index of stock)
                val sizeValue = broadcastSize.value()
                val groupSizeValue = broadcastGroupSize.value()

                // key is key (item in table), value is index in list
                val keyValues = ArrayList<Tuple2<Int, Int>>()
                while (it.hasNext()) {
                    val listIndex = it.next()

                    // for each dimension loop over the other dimensions using addTuples
                    for (dimension in 0 until groupSizeValue) {
                        addTuples(
                            groupSize = groupSizeValue,
                            keyValues = keyValues,
                            value = listIndex,
                            listSize = sizeValue,
                            currentDimension = 0,
                            indexTuple = IntArray(groupSizeValue),
                            skipDimension = dimension,
                        )
                    }
                }
                keyValues.iterator()
            }
        val allPossibleGroups = keys
            .aggregateByKey(
                // each number in table occurs for each dimension as key. The values of those two will be a tuple of (key, stock indices as list)
                // create a new int array to fill up
                zeroValue = getFilledIntArray(
                    size = groupSize,
                    fillValue = -1
                ),

                seqFunction = { base: IntArray, listIndex: Int ->
                    // put listIndex in the first empty spot in base
                    for (i in base.indices) {
                        if (base[i] < 0) {
                            base[i] = listIndex
                            break
                        }
                    }
                    base
                },

                // how to merge partially filled up int arrays
                combineFunction = { a: IntArray, b: IntArray ->
                    // merge a and b
                    var j = 0
                    for (i in a.indices) {
                        if (a[i] < 0) {
                            while (b[j] < 0) {
                                j++
                                if (j == b.size) return@aggregateByKey a
                            }
                            a[i] = b[j]
                            j++
                        }
                    }
                    a
                },
            )
            .values() // finally just take the values

        return createDataset(allPossibleGroups.rdd(), encoder())
    }
}