# Efficiently discovering ternary lagged correlations
Code base used for Thesis: Efficiently discovering ternary lagged correlations, by Jolan J.R. Rensen, TU/e Eindhoven, 2021

The code base is split up into [Main](Main), which houses the STAMP-Pearson-3TS algorithms and no dependencies on Scala or Spark,
and [Spark](Spark) which contains the group-by-key algorithm in context of Spark.
It also contains a compiled version of [ScalaTuplesInKotlin](https://github.com/Jolanrensen/ScalaTuplesInKotlin) in the [lib](lib) folder.

## Examples

For examples regarding the STAMP-Pearson algorithms in Kotlin, take a look [here](Main/src/main/kotlin/nl/jolanrensen/efficientlyDiscoveringTernaryLaggedCorrelations/example/Example.kt). For Java, [here](Main/src/main/java/nl/jolanrensen/efficientlyDiscoveringTernaryLaggedCorrelations/example/Example.java).

For an example of the group-by-key method and how to use STAMP-Pearson-3TS in conjunction with Spark, look [here](Spark/src/main/kotlin/nl/jolanrensen/efficientlyDiscoveringTernaryLaggedCorrelations/example).

## Note

The project is best built using Jetbrains' JDK 11, Gradle 7.2, and needs Kotlin 1.5.32 (due to the Kotlin Spark API's latest version being 1.0.2 atm). 
This also limits the Scala version to 2.12.14 and Spark version to 3.0.0.
