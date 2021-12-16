package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.dataset

import nl.jolanrensen.scalaTuplesInKotlin.first
import nl.jolanrensen.scalaTuplesInKotlin.u
import nl.jolanrensen.thesis.mapInPartitions
import org.apache.spark.sql.Dataset
import org.jetbrains.kotlinx.spark.api.SparkSession
import org.jetbrains.kotlinx.spark.api.col
import org.jetbrains.kotlinx.spark.api.to
import scala.Tuple2
import java.io.Serializable
import java.sql.Timestamp
import kotlin.random.Random

data class StockData(
    val stockExchange: String,
    val stockName: String,
    val instantLongs: LongArray,
    val closingPrices: DoubleArray,
) : Serializable {

    companion object {
        fun read(stocksFileLocation: String, spark: SparkSession, randomize: Boolean, limit: Int?): Dataset<StockData> =
            spark
                .read()
                .parquet(stocksFileLocation)
                .to<StockData>()
                .let {
                    if (randomize)
                        it.mapInPartitions { it u Random.nextInt() }
                            .orderBy(col(Tuple2<*, *>::_2))
                            .mapInPartitions { it.first }
                    else it
                }
                .let {
                    if (limit != null) {
                        it.limit(limit)
                    } else it
                }
    }

    val timeStamps: List<Timestamp>
        get() = instantLongs.map { Timestamp(it) }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as StockData

        if (stockExchange != other.stockExchange) return false
        if (stockName != other.stockName) return false
        if (!closingPrices.contentEquals(other.closingPrices)) return false

        return true
    }

    override fun hashCode(): Int {
        var result = stockExchange.hashCode()
        result = 31 * result + stockName.hashCode()
        result = 31 * result + closingPrices.contentHashCode()
        return result
    }
}