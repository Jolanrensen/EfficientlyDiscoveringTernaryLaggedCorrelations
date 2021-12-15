package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.naive

import org.jetbrains.bio.viktor.F64FlatArray
import java.io.Serializable
import kotlin.math.pow
import kotlin.math.sqrt

object Pearson : AbstractCorrelation(), Serializable {

    override val needsMean = true
    override val needsMin = false
    override val needsMax = false
    override val needsEntropy = false
    override val maxAmountOfSeries = 2


    override fun correlationFunction(
        firstSeries: F64FlatArray,
        secondSeries: F64FlatArray,
        extraArgs: Array<CorrelationExtraData>?
    ): Double {
        val n = firstSeries.length
        require(secondSeries.length == n)

        val x = firstSeries
        val y = secondSeries

        val xMean = extraArgs?.get(0)?.mean ?: x.mean()
        val yMean = extraArgs?.get(1)?.mean ?: y.mean()

        val a = run {
            var sum = 0.0
            var i = 0
            while (i < n) {
                sum += (x[i] - xMean) * (y[i] - yMean)
                i++
            }
            sum
        }

        val b = run {
            var sum = 0.0
            var i = 0
            while (i < n) {
                sum += (x[i] - xMean).pow(2)
                i++
            }
            sum
        }

        val c = run {
            var sum = 0.0
            var i = 0
            while (i < n) {
                sum += (y[i] - yMean).pow(2)
                i++
            }
            sum
        }

        return a / (sqrt(b) * sqrt(c))
    }

    override fun correlationFunction(series: Array<F64FlatArray>, extraArgs: Array<CorrelationExtraData>?): Double {
        val n = series[0].length
        require(series[1].length == n)

        val x = series[0]
        val y = series[1]

        val xMean = extraArgs?.get(0)?.mean ?: x.mean()
        val yMean = extraArgs?.get(1)?.mean ?: y.mean()

        val a = run {
            var sum = 0.0
            var i = 0
            while (i < n) {
                sum += (x[i] - xMean) * (y[i] - yMean)
                i++
            }
            sum
        }

        val b = run {
            var sum = 0.0
            var i = 0
            while (i < n) {
                sum += (x[i] - xMean).pow(2)
                i++
            }
            sum
        }

        val c = run {
            var sum = 0.0
            var i = 0
            while (i < n) {
                sum += (y[i] - yMean).pow(2)
                i++
            }
            sum
        }

        return a / (sqrt(b) * sqrt(c))
    }
}