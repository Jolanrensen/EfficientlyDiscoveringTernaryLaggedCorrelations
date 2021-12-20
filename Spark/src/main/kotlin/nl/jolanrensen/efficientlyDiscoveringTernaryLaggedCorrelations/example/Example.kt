package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.example

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.GroupCalculation.getAllPossibleGroups
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.AggregationWithReducerWithArgWithCompare
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.dataset.StockDataWithoutTime
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3ts
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3ts.AggregationMethod
import nl.jolanrensen.thesis.iteratorOf
import nl.jolanrensen.thesis.mapInPartitions
import nl.jolanrensen.thesis.plusAssign
import org.apache.spark.api.java.JavaSparkContext
import org.apache.spark.broadcast.Broadcast
import org.apache.spark.sql.Dataset
import org.apache.spark.util.BoundedPriorityQueue
import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor.asF64Array
import org.jetbrains.kotlinx.spark.api.*
import java.io.Serializable
import java.time.LocalDate

/**
 * This example shows how to use STAMP-Pearson-3TS correctly in combination with
 * group formation in Apache Spark.
 * (An `object` in Kotlin can be viewed as a singleton)
 */
object Example : Serializable {

    private lateinit var sc: JavaSparkContext

    /** This can store extra data regarding the stocks to save calculations. */
    data class StockExtraInfo(
        val stockExchange: String,
        val stockName: String,
        val closingPrices: DoubleArray,
        val slidingMeans: DoubleArray,
        val slidingStds: DoubleArray,
    ) : Serializable

    /** Index and aggregation method used for the [Result]. */
    data class IndexAgg(val aggregationMethod: String, val a: Int, val b: Int, val c: Int) : Serializable

    /** Data class to store the results in */
    data class Result(
        val stocks: List<StockExtraInfo>,
        val index: IndexAgg,
        val correlation: Double,
    ) : Serializable

    /** Helper function to create a new priority queue for the results */
    @Suppress("NOTHING_TO_INLINE")
    private inline fun newQueue(
        topX: Int,
        reducer: AggregationWithReducerWithArgWithCompare,
    ): BoundedPriorityQueue<Result> =
        BoundedPriorityQueue<Result>(topX) { x, y ->
            reducer.comparator.compare(x.correlation, y.correlation)
        }

    @JvmStatic
    fun main(args: Array<String>) = withSpark(
        master = "local[*]",
        appName = "Example",
    ) {
        sc = JavaSparkContext(spark.sparkContext)

        val windowSize = 10
        val lagBound = 100

        // We've provided some test data used in the thesis
        // This includes pre-processed stocks of 2018
        // Let's read 10 timeseries of it using some helper functions
        val stocksFileLocation = "./data/dataset1.parquet"
        val datesFileLocation = "./data/dataset1dates.parquet"

        val stocksList: List<StockExtraInfo> = StockDataWithoutTime.read(
            stocksFileLocation = stocksFileLocation,
            spark = spark,
            randomize = false,
            limit = 10,
        )
            .mapInPartitions {
                // We will add extra data to save calculations later on
                val (means: F64FlatArray, stds: F64FlatArray) = StampPearson.computeSlidingMeanStd(
                    T = it.closingPrices.asF64Array(),
                    windowSize = windowSize,
                )
                StockExtraInfo(
                    stockExchange = it.stockExchange,
                    stockName = it.stockName,
                    closingPrices = it.closingPrices,
                    slidingMeans = means.data,
                    slidingStds = stds.data,
                )
            }
            .collectAsList()

        val dates: Dataset<LocalDate> = spark
            .read()
            .parquet(datesFileLocation)
            .to<LocalDate>()
            .sort("value")

        val dateList: List<LocalDate> = dates.collectAsList()


        // Since we will distribute the calculations, we are going to broadcast the list
        // to all possible workers
        val stocksListBroadcast: Broadcast<List<StockExtraInfo>> = spark.broadcast(stocksList)

        // Let's generate all possible groups of size 3 for all 30 time series.
        // For this, we use the efficient group-by-key method, described in the thesis
        val groups: Dataset<IntArray> = getAllPossibleGroups(listSize = stocksList.size, groupSize = 3)

        // We'll go over all the groups to get the best aggregation method and window positions for the best
        // correlation value within that group. Then we find the overall top 10 best groups to report
        val results: BoundedPriorityQueue<Result> = groups.mapPartitions { it: Iterator<IntArray> ->
            // we are partitioning the data for each worker and search for the best correlations within groups of 3 time series now
            val stocksList: List<StockExtraInfo> = stocksListBroadcast.value

            // we need an instance of StampPearson per worker
            val stampPearson3ts = StampPearson3ts(maxArraySize = 252) // (dataset1 has series of size 252)

            // Use priority queue to efficiently collect the top 5
            val topXQueue: BoundedPriorityQueue<Result> = newQueue(5, MAX)

            it.forEach { group: IntArray ->
                // the group represents indices for the stocks, so let's get them
                val (tsA: StockExtraInfo, tsB: StockExtraInfo, tsC: StockExtraInfo) = group.map { stocksList[it] }

                // Let's call it!

                var index: IndexAgg? = null
                val bestCorrelation: Double = stampPearson3ts.stampPearson3ts(
                    timeSeriesA = tsA.closingPrices.asF64Array().copy(),
                    timeSeriesB = tsB.closingPrices.asF64Array().copy(),
                    timeSeriesC = tsC.closingPrices.asF64Array().copy(),
                    windowSize = windowSize,
                    reducer = MAX,

                    timeSeriesASlidingMeans = tsA.slidingMeans.asF64Array().copy(),
                    timeSeriesASlidingStds = tsA.slidingStds.asF64Array().copy(),
                    timeSeriesBSlidingMeans = tsB.slidingMeans.asF64Array().copy(),
                    timeSeriesBSlidingStds = tsB.slidingStds.asF64Array().copy(),
                    timeSeriesCSlidingMeans = tsC.slidingMeans.asF64Array().copy(),
                    timeSeriesCSlidingStds = tsC.slidingStds.asF64Array().copy(),

                    lagBound = lagBound,
                ) { aggregation: AggregationMethod, aIndex: Int, bIndex: Int, cIndex: Int ->
                    index = IndexAgg(aggregation.textual, aIndex, bIndex, cIndex)
                }

                // We'll add the result including the index to the queue
                topXQueue.plusAssign(
                    Result(
                        stocks = listOf(tsA, tsB, tsC),
                        index = index!!,
                        correlation = bestCorrelation,
                    )
                )
            }

            return@mapPartitions iteratorOf(topXQueue)
        }
            .reduceK { a: BoundedPriorityQueue<Result>, b: BoundedPriorityQueue<Result> ->
                // let's reduce the queues to get the overall top 10

                val topXQueue: BoundedPriorityQueue<Result> = newQueue(5, MAX)
                topXQueue.plusAssign(a)
                topXQueue.plusAssign(b)

                return@reduceK topXQueue
            }

        // Let's print the results!
        for ((res: List<StockExtraInfo>, resIndex: IndexAgg, corr: Double) in results) {
            val (stock1: StockExtraInfo, stock2: StockExtraInfo, stock3: StockExtraInfo) = res
            val (agg, aIndex, bIndex, cIndex) = resIndex

            println("""
                |Highest correlation $corr found at: 
                |a: ${stock1.stockExchange}, ${stock1.stockName}: ${
                listOf(dateList[aIndex], dateList[aIndex + windowSize - 1])
            } ${stock1.closingPrices.slice(aIndex until aIndex + windowSize)}
                |b: ${stock2.stockExchange}, ${stock2.stockName}: ${
                listOf(dateList[bIndex], dateList[bIndex + windowSize - 1])
            } ${stock2.closingPrices.slice(bIndex until bIndex + windowSize)}
                |c: ${stock3.stockExchange}, ${stock3.stockName}: ${
                listOf(dateList[cIndex], dateList[cIndex + windowSize - 1])
            } ${stock3.closingPrices.slice(cIndex until cIndex + windowSize)}
                |using $agg
            """.trimMargin())

        }


    }


}
