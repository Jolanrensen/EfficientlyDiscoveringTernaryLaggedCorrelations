package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.example;

import kotlin.Pair;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson;
import org.jetbrains.bio.viktor.F64FlatArray;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.F64ArrayHelpersKt.randomF64FlatArray;
import static nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.F64ArrayHelpersKt.toList;

public class Example {

    public Example() {
    }

    public static void main(String[] args) {
        var example = new Example();
        example.matrixProfileExample();
    }


    // TODO

    /**
     * Gives example of how to create a Matrix Profile using this codebase.
     */
    void matrixProfileExample() {

        // let's generate a random time series for the sample
        F64FlatArray timeSeries = randomF64FlatArray(100);

        // first we need an instance of [StampPearson]
        StampPearson stampPearson = new StampPearson(100);

        int windowSize = 10;

        // we'll call it using MAX correlation (so min Euclidean distance)
        Pair<F64FlatArray, int[]> mpResult = stampPearson.matrixProfilePearson(
                timeSeries,
                null,
                windowSize,
                MAX.INSTANCE);

        F64FlatArray matrixProfilePearson = mpResult.getFirst();
        int[] matrixProfileIndices = mpResult.getSecond();

        // if we want the Euclidean distance, we can simply convert the [matrixProfilePearson]
        F64FlatArray matrixProfile = StampPearson.Companion.pearsonCorrelationsToEuclidean(
                matrixProfilePearson,
                windowSize);

        // and we're done!
        System.out.println(
                    "Using time series: "+ toList(timeSeries) +
                    "We get the Matrix Profile: " + toList(matrixProfile) +
                    "With indices: " + List.of(matrixProfileIndices));
    }
}
