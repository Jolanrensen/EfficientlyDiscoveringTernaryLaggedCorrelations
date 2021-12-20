package nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.example;

import kotlin.Pair;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.WindowSkipping;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.aggregations.MAX;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.simplification.RamerDouglasPeucker;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3ts;
import nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.stampPearson.StampPearson3tsWithSkipping;
import org.jetbrains.bio.viktor.F64Array;
import org.jetbrains.bio.viktor.F64FlatArray;

import java.util.*;
import java.util.stream.Collectors;

import static nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.F64ArrayHelpersKt.randomF64FlatArray;
import static nl.jolanrensen.efficientlyDiscoveringTernaryLaggedCorrelations.F64ArrayHelpersKt.toList;
import static org.jetbrains.bio.viktor.F64ArrayKt.asF64Array;

public class Example {

    public Example() {
    }

    public static void main(String[] args) {
        var example = new Example();
        example.matrixProfileExample();
    }

    /**
     * Gives an example for STAMP-Pearson-3TS with 3 time series.
     */
    void stampPearson3tsExample() {

        // let's generate three random time series for this sample
        F64FlatArray tsA = randomF64FlatArray(100);
        F64FlatArray tsB = randomF64FlatArray(100);

        // We can also simply use a double array
        double[] tsCDouble = new double[100];
        {
            var random = new Random();
            for (int i = 0; i < 100; i++) {
                tsCDouble[i] = random.nextDouble();
            }
        }
        F64FlatArray tsC = asF64Array(tsCDouble);

        // we need an instance of [StampPearson3ts]
        StampPearson3ts stampPearson3ts = new StampPearson3ts(100);

        int windowSize = 10;

        // lag bound limits the distance between the windows
        int lagBound = 10;

        // say the sliding means and standard deviations of one of the time series are precalculated somewhere,
        // then we can provide those to the function as well, to save calculations. If not provided, the function will simply
        // calculate them itself.
        var timeSeriesBSlidingMeansStds = StampPearson.Companion.computeSlidingMeanStd(
                tsB,
                windowSize);

        F64FlatArray timeSeriesBSlidingMeans = timeSeriesBSlidingMeansStds.getFirst();
        F64FlatArray timeSeriesBSlidingStds = timeSeriesBSlidingMeansStds.getSecond();

        // if we want to debug and get all the correlation calculations, we can provide a 3D array to the function
        // this does make the calculation slower
        F64Array resultArray = F64Array.Companion.invoke(
                tsC.getLength() - windowSize + 1,
                tsB.getLength() - windowSize + 1,
                tsA.getLength() - windowSize + 1);

        // let's call it!
        double bestCorrelation = stampPearson3ts.stampPearson3ts(
                tsA,
                tsB,
                tsC,
                windowSize,
                MAX.INSTANCE,
                null,
                null,
                timeSeriesBSlidingMeans,
                timeSeriesBSlidingStds,
                null,
                null,
                lagBound,
                resultArray,
                (aggregation, aIndex, bIndex, cIndex) -> {

                    // we can provide a callback as last argument, which will give the indices and aggregation method of the best window positions found

                    System.out.println(
                            "Best index in A is found at: " + aIndex + "\n" +
                                    "Best index in B is found at: " + bIndex + "\n" +
                                    "Best index in C is found at: " + cIndex + "\n" +
                                    "Aggregation method used is: " + aggregation.getTextual()
                    );

                    return null;
                });

        System.out.println("best correlation found is $bestCorrelation");

        System.out.println("all other correlation results are: \n" + resultArray);
    }


    /** Gives an example for STAMP-Pearson-3TS with 3 time series, simplification, and skipping. */
    void stamp3tsWithSkippingExample() {

        // let's generate three random time series for this sample
        List<F64FlatArray> timeSeries = new ArrayList<>() {{
            for (int i = 0; i < 3; i++) {
                add(randomF64FlatArray(100));
            }
        }};

        int windowSize = 10;

        // let's apply simplification to the time series using RamerDouglasPeucker
        double epsilon = 0.01;
        var rdp = new RamerDouglasPeucker(epsilon);
        for (var it : timeSeries) {
            rdp.filterKeepingLengthInPlace(it);
        }

        // Next, let's apply straight line skipping to a copy of the time series, while also keeping the original
        // We must also sort them such that tsA has the fewest NaN values and tsC the most

        class TimeSeries {
            public final F64FlatArray original;
            public final F64FlatArray withSkipping;

            public TimeSeries(F64FlatArray original, F64FlatArray withSkipping) {
                this.original = original;
                this.withSkipping = withSkipping;
            }
        }

        List<TimeSeries> sortedTimeSeries = timeSeries.stream()
                .map(it -> {
                            F64FlatArray withSkipping = it.copy();
                            WindowSkipping.INSTANCE.removeStraightLineContinuity(
                                    withSkipping,
                                    windowSize);

                            return new TimeSeries(it, withSkipping);
                        }
                )
                .sorted(Comparator.comparingInt(a ->
                                WindowSkipping.INSTANCE.numberOfWindowPlacementsSaved(
                                        a.withSkipping,
                                        windowSize)
                        )
                )
                .collect(Collectors.toList());

        TimeSeries tsA = sortedTimeSeries.get(0);
        TimeSeries tsB = sortedTimeSeries.get(1);
        TimeSeries tsC = sortedTimeSeries.get(2);


        // we need an instance of [StampPearson3tsWithSkipping]
        var stampPearson3tsWithSkipping = new StampPearson3tsWithSkipping(
                StampPearson3tsWithSkipping.ContinuitySkippingType.STRAIGHT_LINE,
                100);


        // lag bound limits the distance between the windows
        int lagBound = 10;

        // if we want to debug and get all the correlation calculations, we can provide a 3D array to the function
        // this does make the calculation slower
        F64Array resultArray = F64Array.Companion.invoke(
                tsC.original.getLength() - windowSize + 1,
                tsB.original.getLength() - windowSize + 1,
                tsA.original.getLength() - windowSize + 1);

        // let's call it!
        double bestCorrelation = stampPearson3tsWithSkipping.stampPearson3ts(
                tsA.withSkipping,
                tsB.withSkipping,
                tsC.withSkipping,
                windowSize,
                MAX.INSTANCE,

                tsA.original,
                tsB.original,
                tsC.original,

                null, null, null, null, null, null,

                lagBound,
                resultArray,
                (aggregation, aIndex, bIndex, cIndex) -> {
                    // we can provide a callback as last argument, which will give the indices and aggregation method of the best window positions found

                    System.out.println("Best index in A is found at: " + aIndex + "\n" +
                            "Best index in B is found at: " + bIndex + "\n" +
                            "Best index in C is found at: " + cIndex + "\n" +
                            "Aggregation method used is: " + aggregation.getTextual());

                    return null;
                });

        System.out.println("best correlation found is " + bestCorrelation);

        System.out.println("all other correlation results are: \n" + resultArray);
    }

    /** Gives example of how to use STAMP-Pearson with 2 time series using this codebase. */
    void stampPearsonExample() {

        // let's generate two random time series for this sample
        F64FlatArray tsA = randomF64FlatArray(100);
        F64FlatArray tsB = randomF64FlatArray(100);

        // first we need an instance of [StampPearson]
        var stampPearson = new StampPearson(100);

        int windowSize = 10;

        // if we want to debug and get all the correlation calculations, we can provide a 2D array to the function
        // this does make the calculation slower
        F64Array resultArray = F64Array.Companion.invoke(
                tsB.getLength() - windowSize + 1,
                tsA.getLength() - windowSize + 1);

        double bestCorrelation = stampPearson.stampPearson(
                tsA.copy(), // it's good practice providing a copy to the array, since they are likely overwritten by the algorithm
                tsB.copy(),
                windowSize,
                MAX.INSTANCE,
                null, null, null, null, null, null,
                resultArray,
                (aIndex, bIndex) -> {
                    // we can provide a callback as last argument, which will give the indices of the best window positions found

                    System.out.println("Best index in A is found at: " + aIndex + "\n" +
                            "Best index in B is found at: " + bIndex);
                    return null;
                });

        System.out.println("best correlation found is " + bestCorrelation);

        System.out.println("all other correlation results are:" +
                Arrays.stream(resultArray.toGenericArray())
                        .map(it -> List.of((double[]) it).toString())
                        .collect(Collectors.joining("\n"))
        );
    }

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
                "Using time series: " + toList(timeSeries) + "\n" +
                        "We get the Matrix Profile: " + toList(matrixProfile) + "\n" +
                        "With indices: " + List.of(matrixProfileIndices));
    }
}
