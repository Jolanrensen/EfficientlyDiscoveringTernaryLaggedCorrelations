package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.naive

data class CorrelationExtraData(
    val mean: Double? = null,
    val min: Double? = null,
    val max: Double? = null,
    val entropy: Double? = null,
)
