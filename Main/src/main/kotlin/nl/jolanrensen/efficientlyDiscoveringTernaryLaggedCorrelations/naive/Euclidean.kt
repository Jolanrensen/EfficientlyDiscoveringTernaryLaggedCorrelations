package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.naive

import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.powInPlace
import org.jetbrains.bio.viktor.F64FlatArray
import java.io.Serializable
import kotlin.math.sqrt

object Euclidean : AbstractCorrelation(), Serializable {

    override val needsMean = false
    override val needsMin = false
    override val needsMax = false
    override val needsEntropy = false
    override val maxAmountOfSeries = 2


    override fun correlationFunction(series: Array<F64FlatArray>, extraArgs: Array<CorrelationExtraData>?): Double {
        checkRequirements(series, extraArgs)
        val (x, y) = series

        return sqrt(
            (x - y)
                .powInPlace(2.0)
                .sum()
        )
    }

    override fun correlationFunction(
        firstSeries: F64FlatArray,
        secondSeries: F64FlatArray,
        extraArgs: Array<CorrelationExtraData>?,
    ): Double {
        checkRequirements(arrayOf(firstSeries, secondSeries), extraArgs)

        return sqrt(
            (firstSeries - secondSeries)
                .powInPlace(2.0)
                .sum()
        )
    }
}