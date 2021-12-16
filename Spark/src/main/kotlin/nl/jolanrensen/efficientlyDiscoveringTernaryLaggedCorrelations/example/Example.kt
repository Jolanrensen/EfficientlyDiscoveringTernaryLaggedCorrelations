package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.example

import org.apache.spark.api.java.JavaSparkContext
import org.jetbrains.kotlinx.spark.api.sparkContext
import org.jetbrains.kotlinx.spark.api.withSpark
import java.io.Serializable

object Example : Serializable {

    private lateinit var sc: JavaSparkContext

    @JvmStatic
    fun main(args: Array<String>) = withSpark(
        master = "local[*]",
        appName = "Example",
    ) {
        sc = JavaSparkContext(spark.sparkContext)

        println("test")


    }


}
