import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    kotlin("jvm") version "1.6.10"
    application
    java
}

group = "nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
    mavenCentral()
    maven(url = "https://jitpack.io")
    maven(url = "https://maven.pkg.jetbrains.space/mipt-npm/p/sci/maven")
    maven(url = "https://maven.pkg.jetbrains.space/data2viz/p/maven/public")
}

val scalaVersion = "2.12.14"
val sparkVersion = "3.0.0"

dependencies {

    implementation(project(":Main"))

    implementation(kotlin("stdlib"))
    testImplementation("org.junit.jupiter:junit-jupiter-api:5.8.2")
    testRuntimeOnly("org.junit.jupiter:junit-jupiter-engine")

    // Scala libraries. Version must match Kotlin Spark API
    implementation(group = "org.scala-lang", name = "scala-library", version = scalaVersion) // <4>
    implementation(group = "org.scala-lang", name = "scala-reflect", version = scalaVersion)

    implementation(
        group = "org.jetbrains.kotlinx.spark",
        name = "kotlin-spark-api-3.0",
        version = "1.0.2",
    )

    implementation(files("../lib/nl/jolanrensen/scalaTuplesInKotlin/1.0-SNAPSHOT/scalaTuplesInKotlin-1.0-SNAPSHOT.jar"))

    // Apache Spark
    implementation(group = "org.apache.spark", name = "spark-sql_2.12", version = sparkVersion)
    implementation(group = "org.apache.spark", name = "spark-core_2.12", version = sparkVersion)
    implementation(group = "org.apache.spark", name = "spark-streaming_2.12", version = sparkVersion)

    // Some NumPy ndarray features, uses SIMD
    implementation(group = "org.jetbrains.bio", name = "viktor", version = "1.2.0")
    implementation(group = "org.slf4j", name = "slf4j-log4j12", version = "1.7.25")

    implementation(group = "org.jetbrains.kotlinx", name = "kotlinx-coroutines-core", version = "1.4.3")
    testImplementation(group = "org.jetbrains.kotlinx", name = "kotlinx-coroutines-test", version = "1.4.3")

    implementation(group = "org.jetbrains.kotlinx", name = "kotlinx-datetime", version = "0.3.0")


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
    languageVersion = "1.6"
}

val compileKotlin: KotlinCompile by tasks
compileKotlin.applyOptions()

val compileTestKotlin: KotlinCompile by tasks
compileTestKotlin.applyOptions()