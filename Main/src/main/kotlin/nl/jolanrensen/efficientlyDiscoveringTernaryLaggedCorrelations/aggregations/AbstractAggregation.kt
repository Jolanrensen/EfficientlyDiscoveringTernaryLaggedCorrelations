package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations

import org.jetbrains.bio.viktor.F64FlatArray
import org.jetbrains.bio.viktor.asF64Array
import java.io.Serializable

abstract class AbstractAggregation : Serializable, (Array<F64FlatArray>) -> F64FlatArray {

    override fun invoke(p1: Array<F64FlatArray>): F64FlatArray = aggregationFunction(p1)

    @JvmName("aggregationFunctionArray")
    fun aggregationFunction(
        seriesArray: Array<F64FlatArray>,
    ): F64FlatArray = aggregationFunction(*seriesArray)

    abstract fun aggregationFunction(
        vararg series: F64FlatArray
    ): F64FlatArray

    @JvmName("aggregationFunctionArray")
    fun aggregationFunction(
        seriesArray: Array<DoubleArray>
    ): DoubleArray = aggregationFunction(*seriesArray)

    fun aggregationFunction(
        vararg seriesArray: DoubleArray
    ): DoubleArray = aggregationFunction(
        Array(seriesArray.size) { seriesArray[it].asF64Array() }
    ).toDoubleArray()

    @Throws(IllegalArgumentException::class)
    protected fun checkRequirements(series: Array<F64FlatArray>) {
        val n = series.first().length
        require(series.all { it.length == n }) { "All series should have equal lengths." }
    }
}