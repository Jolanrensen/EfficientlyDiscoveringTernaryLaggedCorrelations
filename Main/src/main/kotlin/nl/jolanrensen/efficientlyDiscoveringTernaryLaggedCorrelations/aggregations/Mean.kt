package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations

import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.F64FlatArray


object Mean : AbstractAggregation() {

    @Suppress("UNCHECKED_CAST")
    override fun aggregationFunction(vararg series: F64FlatArray): F64FlatArray {
        checkRequirements(series as Array<F64FlatArray>)

        return F64Array(series.size, series.first().length) { i, j -> series[i][j] }
            .along(1)
            .map { it.mean() }
            .let { seq ->
                val iter = seq.iterator()
                F64Array(seq.count()) { iter.next() }
            } as F64FlatArray
    }
}