package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.naive

import org.jetbrains.bio.viktor.F64FlatArray
import java.io.Serializable

abstract class AbstractCorrelation : Serializable {

    abstract fun correlationFunction(
        series: Array<F64FlatArray>,
        extraArgs: Array<CorrelationExtraData>? = null,
    ): Double

    abstract fun correlationFunction(
        firstSeries: F64FlatArray,
        secondSeries: F64FlatArray,
        extraArgs: Array<CorrelationExtraData>? = null,
    ): Double

    fun correlationFunction(
        firstSeries: F64FlatArray,
        vararg otherSeries: F64FlatArray,
        extraArgs: Array<CorrelationExtraData>? = null,
    ): Double = correlationFunction(
        series = arrayOf(firstSeries, *otherSeries),
        extraArgs = extraArgs,
    )

    // Return true of the correlation needs these values in the CorrelationExtraData
    abstract val needsMean: Boolean
    abstract val needsMin: Boolean
    abstract val needsMax: Boolean
    abstract val needsEntropy: Boolean

    // Null if there isn't a max
    abstract val maxAmountOfSeries: Int?

    @Throws(IllegalArgumentException::class)
    protected fun checkRequirements(series: Array<F64FlatArray>, extraArgs: Array<CorrelationExtraData>?) {
        val n = series.first().length
        require(series.all { it.length == n }) { "All series should have equal lengths." }

        maxAmountOfSeries?.let {
            require(series.size <= it) { "Max amount of time series for this correlation function is $maxAmountOfSeries." }
        }

        if (needsMean) require(extraArgs!!.all { it.mean != null })
        if (needsMin) require(extraArgs!!.all { it.min != null })
        if (needsMax) require(extraArgs!!.all { it.max != null })
        if (needsEntropy) require(extraArgs!!.all { it.entropy != null })
    }
}