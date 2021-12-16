package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.simplification


import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.f64FlatArrayOf
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.indices
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.invoke
import org.jetbrains.bio.viktor.F64FlatArray
import space.kscience.kmath.functions.asFunction
import space.kscience.kmath.interpolation.LinearInterpolator
import space.kscience.kmath.interpolation.interpolatePolynomials
import space.kscience.kmath.operations.RealField
import space.kscience.kmath.structures.asBuffer
import java.io.Serializable
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt

internal class RamerDouglasPeucker(
    epsilon: Double,
) : Serializable {

    init {
        require(epsilon >= 0) { "Epsilon must be > 0" }
    }

    var epsilon: Double = epsilon
        set(value) {
            require(epsilon >= 0) { "Epsilon must be > 0" }
            field = value
        }


    fun filter(data: F64FlatArray): Pair<F64FlatArray, IntArray> = ramerDouglasPeuckerFunction(
        points = data,
        indices = data.indices.toList().toIntArray(),
        startIndex = 0,
        endIndex = data.length - 1,
    )

    private fun ramerDouglasPeuckerFunction(
        points: F64FlatArray,
        indices: IntArray,
        startIndex: Int,
        endIndex: Int,
    ): Pair<F64FlatArray, IntArray> {
        var dmax = 0.0
        var idx = 0
        val a = (endIndex - startIndex).toDouble()
        val b = points[endIndex] - points[startIndex]
        val c = -(b * startIndex - a * points[startIndex])
        val norm = sqrt(a.pow(2.0) + b.pow(2.0))
        for (i in startIndex + 1 until endIndex) {
            val distance = abs(b * i - a * points[i] + c) / norm
            if (distance > dmax) {
                idx = i
                dmax = distance
            }
        }
        if (dmax >= epsilon) {
            val (recursiveResult1, recursiveResult1Indices) = ramerDouglasPeuckerFunction(
                points = points,
                indices = indices,
                startIndex = startIndex,
                endIndex = idx,
            )
            val (recursiveResult2, recursiveResult2Indices) = ramerDouglasPeuckerFunction(
                points = points,
                indices = indices,
                startIndex = idx,
                endIndex = endIndex,
            )
            val result = F64FlatArray(
                recursiveResult1.length - 1 + recursiveResult2.length
            )
            val resultIndices = IntArray(
                recursiveResult1Indices.size - 1 + recursiveResult2Indices.size
            )

            System.arraycopy(
                recursiveResult1.data,
                0,
                result.data,
                0,
                recursiveResult1.length - 1
            )
            System.arraycopy(
                recursiveResult2.data,
                0,
                result.data,
                recursiveResult1.length - 1,
                recursiveResult2.length
            )

            System.arraycopy(
                recursiveResult1Indices,
                0,
                resultIndices,
                0,
                recursiveResult1Indices.size - 1
            )
            System.arraycopy(
                recursiveResult2Indices,
                0,
                resultIndices,
                recursiveResult1Indices.size - 1,
                recursiveResult2Indices.size
            )

            return result to resultIndices
        } else {

            return f64FlatArrayOf(points[startIndex], points[endIndex]) to
                    intArrayOf(indices[startIndex], indices[endIndex])
        }
    }

    fun filterKeepingLengthInPlace(data: F64FlatArray) {
        if (epsilon == 0.0) return

        val (fPoints, fIndices) = filter(data)
        val interpolator = LinearInterpolator(RealField).interpolatePolynomials(
            x = fIndices.map { it.toDouble() }.asBuffer(),
            y = fPoints.toDoubleArray().asBuffer(),
        ).asFunction(RealField)

        for (i in data.indices) {
            data[i] = interpolator(i.toDouble()) ?: data[i]
        }
    }

    fun filterKeepingLength(data: F64FlatArray): F64FlatArray {
        val result = data.copy()
        filterKeepingLengthInPlace(result)
        return result
    }

}