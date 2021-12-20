import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    kotlin("jvm") version "1.5.32"
    application
    java
}

group = "nl.jolanrensen.effecientlyDiscoveryingTernaryLaggedCorrelations"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
    mavenCentral()
    maven(url = "https://jitpack.io")
    maven(url = "https://maven.pkg.jetbrains.space/mipt-npm/p/sci/maven")
    maven(url = "https://maven.pkg.jetbrains.space/data2viz/p/maven/public")
}

dependencies {
    implementation(kotlin("stdlib"))
    testImplementation("org.junit.jupiter:junit-jupiter-api:5.8.2")
    testRuntimeOnly("org.junit.jupiter:junit-jupiter-engine")

    // Some NumPy ndarray features, uses SIMD
    implementation(group = "org.jetbrains.bio", name = "viktor", version = "1.2.0")
    implementation(group = "org.slf4j", name = "slf4j-log4j12", version = "1.7.25")

    implementation(group = "org.jetbrains.kotlinx", name = "kotlinx-coroutines-core", version = "1.4.3")
    testImplementation(group = "org.jetbrains.kotlinx", name = "kotlinx-coroutines-test", version = "1.4.3")

    // FFT https://sites.google.com/site/piotrwendykier/software/jtransforms
    implementation(
        group = "com.github.wendykierp",
        name = "JTransforms",
        version = "3.1",
    )

    // regression etc
    implementation(group = "org.nield", name = "kotlin-statistics", version = "1.2.1")

    // Arrays, matrices etc
    val kmathGroup = "space.kscience"
    val kMathVersion = "0.2.1"
    api(group = kmathGroup, name = "kmath-core", version = kMathVersion)
    api(group = kmathGroup, name = "kmath-commons", version = kMathVersion)
    api(group = kmathGroup, name = "kmath-complex", version = kMathVersion)
    api(group = kmathGroup, name = "kmath-for-real", version = kMathVersion)
    api(group = kmathGroup, name = "kmath-dimensions", version = kMathVersion)
    api(group = kmathGroup, name = "kmath-stat", version = kMathVersion)

}

tasks.getByName<Test>("test") {
    useJUnitPlatform()
}

fun KotlinCompile.applyOptions() = kotlinOptions {
    freeCompilerArgs = listOf(
        "-Xinline-classes",
        "-Xopt-in=org.mylibrary.OptInAnnotation",
        "-Xopt-in=org.jetbrains.numkt.core.ExperimentalNumkt",
        "-Xopt-in=kotlin.RequiresOptIn",
        "-Xunrestricted-builder-inference",
        "-Xopt-in=kotlin.time.ExperimentalTime", // measureTimed etc
        "-Xopt-in=kotlin.ExperimentalStdlibApi",  // buildList etc
    )
    jvmTarget = "11"
    languageVersion = "1.5"
}

val compileKotlin: KotlinCompile by tasks
compileKotlin.applyOptions()

val compileTestKotlin: KotlinCompile by tasks
compileTestKotlin.applyOptions()