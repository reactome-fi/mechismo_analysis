/*
 * Created on Apr 6, 2007
 *
 */
package org.reactome.r3.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.HypergeometricDistribution;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math.stat.inference.TestUtils;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Assert;
import org.junit.Test;

import cern.jet.random.Binomial;
import cern.jet.random.ChiSquare;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;

public class MathUtilities {

    //private static DistributionFactory distFactory;
    private static RandomEngine randomeEngine;

    static {
        //distFactory = DistributionFactory.newInstance();
        randomeEngine = new DRand();
    }

    //from http://bioconductor.org/packages/release/bioc/vignettes/EmpiricalBrownsMethod/inst/doc/ebmVignette.html
    private final Double[] glypPvals = {1.990673e-01, 3.111716e-10, 9.354057e-01, 7.620790e-01, 8.072350e-03, 7.798301e-02, 2.564484e-01};
    private final Double[] row1 = {6.126, 6.919, 6.099, 7.381, 5.565, 5.732, 7.851, 6.47, 6.073, 6.117, 6.114, 6.115, 6.403, 6.327, 7.162, 6.87, 7.2, 6.321, 5.942, 5.811, 4.241, 4.418, 4.893, 7.782, 6.426, 6.579, 6.134, 7.342, 6.956, 6.292, 6.377, 5.995, 6.219, 5.013, 5.738, 4.725, 5.54, 8.766, 4.741, 6.327, 4.493, 4.333, 6.201, 5.028, 3.907, 6.376, 8.309, 9.41, 8.17, 6.718, 6.516, 5.453, 6.458, 6.389, 4.57, 6.338, 6.556, 5.488, 5.559, 6.707, 5.499, 5.16, 6.298, 6.336, 6.972, 6.328, 5.703, 6.274, 5.535, 5.419, 4.553, 4.25, 5.35, 6.406, 4.145, 5.886, 7.227, 7.235, 10.0, 6.963, 7.408, 7.29, 5.722, 5.769, 4.826, 7.041, 7.184, 5.998, 5.817, 6.082, 5.41, 6.929, 4.515, 7.006, 7.112, 4.422, 6.408, 6.856, 7.359, 6.224, 6.905, 6.486, 5.484, 4.3, 6.518, 5.344, 7.599, 6.432, 5.34, 5.226, 5.859, 6.195, 5.003, 5.403, 4.685, 7.215, 6.995, 6.384, 3.65, 4.731, 5.378, 4.174, 4.123, 7.116, 5.759, 5.376, 6.901, 8.256, 7.716, 8.622, 5.649, 6.079, 6.015, 6.501, 7.094, 4.574, 5.617, 6.252, 6.1, 6.5, 5.052, 6.717, 5.621, 6.679, 4.986, 8.053, 5.275, 4.203, 5.085, 5.582, 5.061, 4.87, 2.579, 5.104, 5.579, 4.983, 6.475, 7.277, 6.014, 6.585, 5.289, 6.361, 5.23, 5.761, 6.042, 6.489, 5.534};
    private final Double[] row2 = {9.575, 7.808, 9.635, 11.797, 7.69, 9.197, 12.216, 6.528, 8.054, 7.882, 7.886, 8.882, 8.72, 9.559, 8.677, 9.214, 9.538, 8.109, 9.651, 8.575, 7.32, 9.709, 6.988, 10.596, 8.871, 7.849, 9.167, 8.906, 9.13, 8.081, 9.378, 8.09, 9.734, 8.534, 8.091, 8.218, 9.202, 8.979, 7.664, 9.536, 8.603, 8.978, 10.002, 7.909, 9.634, 9.696, 8.601, 10.072, 12.189, 9.71, 8.285, 8.564, 7.458, 8.888, 9.083, 9.91, 9.325, 9.645, 8.218, 9.109, 9.008, 12.338, 10.663, 8.597, 10.417, 8.422, 6.875, 9.483, 8.186, 7.761, 8.119, 7.275, 10.579, 9.118, 9.947, 9.505, 8.66, 7.577, 13.355, 9.059, 7.94, 9.998, 6.515, 8.684, 9.05, 9.57, 9.397, 8.044, 10.918, 13.106, 11.003, 9.976, 8.075, 6.977, 11.712, 8.398, 10.304, 10.052, 8.597, 9.512, 10.039, 8.378, 10.002, 8.382, 9.212, 6.759, 8.682, 9.62, 8.318, 8.764, 8.73, 10.299, 8.322, 8.656, 7.706, 9.543, 9.011, 9.057, 7.0, 9.386, 7.975, 8.909, 9.463, 9.295, 8.45, 9.937, 6.65, 10.127, 9.746, 8.296, 7.311, 9.318, 8.934, 8.805, 9.789, 9.413, 8.37, 10.21, 9.141, 6.955, 8.079, 9.589, 6.831, 9.085, 9.154, 11.655, 9.082, 9.223, 8.597, 8.989, 10.434, 9.642, 9.232, 8.814, 8.434, 9.578, 8.82, 8.07, 8.322, 8.746, 10.935, 6.726, 10.052, 7.891, 9.328, 7.602, 7.856};
    private final Double[] row3 = {2.405, 1.93, 0.502, 0.0, 0.689, 5.311, 0.009, 0.607, 1.139, 0.013, 4.765, 2.273, 0.013, 3.916, 0.003, 1.332, 0.009, 0.579, 0.583, 2.281, 3.281, 4.774, 1.207, 1.112, 5.355, 0.011, 1.321, 0.011, 0.004, 2.388, 3.911, 1.182, 6.815, 1.708, 0.0, 0.01, 2.509, 0.006, 0.001, 3.012, 1.039, 1.934, 9.058, 1.746, 0.01, 4.903, 0.01, 0.007, 0.579, 2.438, 2.728, 3.351, 1.64, 0.549, 6.566, 6.052, 0.004, 3.398, 0.589, 3.474, 1.242, 0.001, 3.406, 0.009, 1.566, 0.013, 0.007, 2.134, 0.898, 1.546, 0.794, 0.009, 8.004, 3.989, 6.051, 5.209, 1.898, 0.646, 0.011, 0.0, 1.323, 2.117, 0.004, 5.216, 1.366, 0.53, 1.951, 0.011, 3.058, 0.001, 4.6, 6.518, 0.007, 0.004, 2.481, 3.457, 3.211, 1.035, 1.966, 5.214, 0.0, 0.006, 0.004, 2.606, 2.752, 0.87, 2.683, 0.003, 0.014, 3.378, 2.676, 2.029, 4.449, 5.444, 2.301, 3.579, 0.577, 0.006, 0.006, 7.774, 3.398, 3.287, 4.486, 0.009, 5.027, 2.444, 1.245, 0.902, 1.629, 0.013, 0.003, 1.567, 5.476, 0.013, 1.45, 2.533, 0.634, 1.337, 4.369, 0.013, 0.006, 3.383, 1.503, 1.866, 4.222, 0.006, 5.397, 0.964, 2.728, 10.351, 1.92, 5.374, 4.9, 1.726, 2.213, 5.718, 2.076, 0.006, 1.3, 3.407, 2.811, 1.088, 2.934, 4.228, 1.396, 1.487, 1.842};
    private final Double[] row4 = {7.654, 7.815, 5.285, 4.633, 4.224, 3.633, 6.361, 3.21, 3.895, 6.641, 5.513, 5.984, 4.626, 3.688, 4.274, 4.504, 8.284, 4.331, 4.332, 4.969, 6.937, 3.943, 5.153, 9.903, 5.663, 5.721, 5.812, 6.445, 8.847, 4.123, 4.315, 5.426, 2.211, 2.742, 6.201, 3.778, 4.467, 3.065, 5.025, 5.164, 4.153, 3.744, 3.896, 4.964, 3.744, 4.301, 4.751, 13.051, 5.026, 5.317, 5.28, 4.279, 5.349, 8.514, 3.583, 6.793, 7.447, 2.75, 4.449, 5.422, 5.279, 5.657, 3.581, 6.852, 7.808, 7.144, 4.416, 4.774, 3.122, 1.836, 1.964, 3.726, 4.109, 5.702, 3.48, 2.871, 3.939, 6.492, 11.632, 3.493, 3.931, 10.311, 6.331, 4.342, 2.971, 5.28, 5.121, 4.253, 1.944, 5.298, 5.308, 3.345, 3.21, 5.526, 6.278, 4.951, 4.524, 4.651, 2.959, 5.773, 9.257, 4.275, 5.484, 5.198, 6.542, 6.153, 7.285, 5.86, 4.149, 3.917, 6.506, 5.589, 3.271, 4.774, 3.237, 6.781, 6.076, 6.143, 2.465, 5.025, 5.225, 4.174, 4.081, 7.344, 4.029, 3.888, 6.194, 9.809, 8.86, 11.037, 4.47, 5.948, 5.441, 4.522, 4.796, 3.545, 5.814, 4.937, 4.596, 6.432, 3.962, 8.346, 4.076, 7.088, 7.152, 3.342, 4.038, 2.389, 5.712, 5.172, 7.197, 2.627, 2.801, 5.22, 3.138, 4.57, 5.953, 7.624, 6.97, 6.302, 4.289, 4.689, 3.454, 6.618, 3.638, 7.48, 2.772};
    private final Double[] row5 = {4.81, 9.692, 3.173, 3.69, 3.343, 3.483, 3.142, 3.564, 5.169, 3.592, 5.545, 3.762, 3.127, 1.865, 3.347, 3.471, 3.354, 3.063, 3.595, 5.203, 4.241, 3.509, 3.245, 4.0, 3.226, 3.867, 3.174, 5.061, 3.84, 3.245, 3.493, 3.589, 1.493, 4.535, 5.029, 3.074, 4.613, 3.934, 1.966, 4.537, 3.31, 3.682, 4.207, 3.567, 4.308, 4.027, 3.539, 5.259, 3.868, 2.845, 4.889, 3.889, 5.33, 6.935, 4.221, 3.849, 2.635, 5.292, 3.078, 3.593, 3.632, 3.492, 3.737, 2.961, 4.282, 2.589, 3.733, 2.956, 4.11, 3.98, 3.479, 3.414, 3.745, 2.69, 3.886, 3.897, 4.597, 5.264, 3.404, 4.387, 4.46, 4.147, 5.171, 5.903, 3.314, 4.144, 4.473, 2.973, 2.909, 2.625, 3.152, 3.881, 2.76, 2.614, 5.469, 3.79, 3.41, 4.525, 2.661, 4.092, 2.35, 3.837, 4.365, 3.159, 2.858, 6.084, 4.218, 5.314, 2.813, 3.009, 4.006, 4.253, 3.489, 3.135, 4.246, 5.589, 3.362, 5.383, 3.761, 4.526, 3.399, 3.536, 2.802, 4.391, 3.115, 2.763, 4.226, 3.844, 4.781, 3.796, 4.645, 3.558, 4.251, 3.681, 3.931, 2.859, 3.524, 4.596, 4.097, 4.034, 3.865, 4.176, 3.428, 4.627, 4.389, 4.4, 3.602, 4.203, 4.504, 2.306, 3.848, 3.011, 1.58, 4.479, 4.183, 2.54, 5.684, 9.363, 3.716, 5.467, 5.251, 5.256, 4.138, 4.651, 3.975, 5.552, 3.931};
    private final Double[] row6 = {10.272, 10.072, 9.396, 9.826, 10.463, 9.476, 9.938, 10.013, 11.087, 10.351, 10.309, 10.256, 9.934, 10.08, 10.085, 10.03, 10.045, 10.04, 9.964, 9.848, 10.1, 9.774, 10.244, 10.436, 9.768, 10.161, 9.863, 10.156, 9.684, 9.944, 9.149, 9.704, 9.616, 10.277, 9.724, 9.914, 9.979, 9.813, 9.819, 10.198, 9.874, 10.055, 10.252, 9.953, 9.821, 9.785, 10.084, 8.93, 9.962, 9.629, 10.166, 10.917, 10.032, 9.865, 9.536, 9.978, 10.555, 9.297, 9.978, 10.139, 10.203, 10.333, 9.75, 9.584, 10.6, 10.274, 9.77, 9.898, 9.346, 9.202, 9.427, 9.836, 9.421, 9.846, 9.67, 9.899, 9.633, 9.677, 9.847, 9.577, 10.642, 10.067, 10.086, 10.136, 9.449, 10.103, 10.266, 9.817, 10.047, 10.05, 9.6, 9.928, 10.024, 9.942, 10.354, 9.877, 10.113, 10.694, 9.857, 9.173, 9.754, 9.384, 9.578, 9.926, 10.193, 9.669, 9.71, 10.131, 9.69, 9.774, 9.829, 9.518, 9.839, 9.704, 9.229, 9.665, 9.965, 9.568, 9.915, 10.108, 9.856, 10.412, 9.523, 9.829, 9.473, 9.608, 10.07, 9.863, 10.243, 9.606, 10.246, 9.674, 11.011, 10.0, 9.881, 9.961, 10.241, 10.363, 9.7, 9.538, 9.828, 10.193, 9.715, 9.953, 9.961, 10.113, 9.687, 10.217, 9.402, 9.752, 9.537, 9.705, 8.98, 11.078, 9.984, 9.804, 10.216, 10.135, 10.015, 9.881, 11.13, 10.37, 9.888, 9.948, 10.304, 10.345, 10.472};
    private final Double[] row7 = {10.587, 12.112, 10.195, 9.667, 10.559, 10.387, 9.543, 11.559, 10.514, 10.213, 11.223, 10.496, 10.009, 10.541, 10.98, 10.06, 9.003, 10.662, 10.551, 11.213, 10.725, 9.962, 9.359, 9.615, 12.226, 11.025, 10.338, 10.339, 10.117, 9.868, 9.852, 10.5, 10.827, 10.259, 10.759, 11.147, 10.474, 10.005, 10.153, 9.879, 9.578, 10.544, 11.689, 10.163, 10.161, 10.121, 10.396, 10.323, 10.449, 10.932, 11.631, 9.971, 11.393, 10.136, 10.545, 10.188, 9.224, 9.27, 10.344, 10.078, 10.483, 10.376, 10.051, 9.885, 8.815, 11.068, 10.245, 9.831, 9.854, 9.673, 10.649, 10.236, 10.478, 11.714, 10.822, 11.476, 9.457, 9.768, 10.013, 10.862, 9.98, 9.911, 8.728, 9.941, 10.389, 9.715, 10.014, 11.087, 10.087, 10.059, 11.094, 10.658, 9.776, 10.729, 9.933, 10.709, 9.886, 10.017, 10.429, 10.64, 9.958, 11.068, 9.806, 10.603, 10.596, 8.971, 10.564, 9.531, 10.167, 10.274, 9.93, 9.984, 10.24, 10.65, 9.867, 10.46, 9.817, 11.022, 10.497, 11.016, 11.296, 10.432, 10.613, 10.15, 10.006, 10.531, 10.266, 10.24, 10.178, 11.717, 9.379, 10.306, 10.743, 10.748, 10.787, 10.033, 10.271, 9.839, 10.75, 10.609, 9.467, 9.916, 10.056, 10.77, 12.154, 9.765, 10.233, 9.97, 10.664, 11.994, 10.446, 10.265, 10.587, 10.692, 10.56, 11.166, 11.562, 11.914, 10.771, 10.779, 10.23, 10.034, 9.158, 9.487, 10.209, 11.78, 12.436};

    /**
     * It seems that there are some bugs in the package for LogRankTest(). The results
     * from this implementation is differnet from R.
     *
     * @return
     */
//    public static double logRankSurvivalTest(List<Double> time1,
//                                             List<Double> censor1,
//                                             List<Double> time2,
//                                             List<Double> censor2) {
//        // Convert an double list to a single array
//        double[] time11 = convertDoubleListToArray(time1);
//        double[] censor11 = convertDoubleListToArray(censor1);
//        double[] time21 = convertDoubleListToArray(time2);
//        double[] censor21 = convertDoubleListToArray(censor2);
//        LogRankTest test = new LogRankTest(time11, censor11, time21, censor21);
//        return test.pValue;
//    }
    private static double[] convertDoubleListToArray(List<Double> list) {
        double[] rtn = new double[list.size()];
        for (int i = 0; i < list.size(); i++)
            rtn[i] = list.get(i);
        return rtn;
    }

    public static double calculateDistance(List<Integer> v1,
                                           List<Integer> v2) {
        if (v1.size() != v2.size())
            throw new IllegalArgumentException("The passed two boolean vectors should have the same length!");
        double t = 0.0d;
        for (int i = 0; i < v1.size(); i++) {
            int i1 = v1.get(i);
            int i2 = v2.get(i);
            t += (double) (i1 - i2) * (i1 - i2);
        }
        return Math.sqrt(t);
    }

    public static double log2(double value) {
        return Math.log(value) / Math.log(Math.E);
    }

    public static double calculateHammingDistance(List<Boolean> vector1,
                                                  List<Boolean> vector2) {
        if (vector1.size() != vector2.size())
            throw new IllegalArgumentException("The passed two boolean vectors should have the same length!");
        int dist = 0;
        for (int i = 0; i < vector1.size(); i++) {
            Boolean b1 = vector1.get(i);
            Boolean b2 = vector2.get(i);
            if (!b1.equals(b2))
                dist++;
        }
        return dist;
    }

    public static double calculateNetworkDistance(List<Boolean> vector1,
                                                  List<Boolean> vector2,
                                                  List<Set<String>> clusters) {
        if (vector1.size() != vector2.size())
            throw new IllegalArgumentException("The passed two boolean vectors should have the same length!");
        int similarity = 0;
        int length = vector1.size();
        int total = 0;
        for (int i = 0; i < vector1.size(); i++) {
            Set<String> cluster = clusters.get(i);
            int weight = cluster.size();
            total += weight;
            Boolean b1 = vector1.get(i);
            Boolean b2 = vector2.get(i);
            if (b1.equals(b2))
                similarity += weight;
        }
        return total - similarity;
    }

    public static double calculateTTest(List<Double> sample1,
                                        List<Double> sample2) throws MathException {
        double[] sampleArray1 = new double[sample1.size()];
        for (int i = 0; i < sample1.size(); i++)
            sampleArray1[i] = sample1.get(i);
        double[] sampleArray2 = new double[sample2.size()];
        for (int i = 0; i < sample2.size(); i++)
            sampleArray2[i] = sample2.get(i);
        return TestUtils.tTest(sampleArray1, sampleArray2);
    }

    public static double calculateMean(List<Double> values) {
        double total = 0.0d;
        for (Double value : values)
            total += value;
        return total / values.size();
    }

    public static double calculateJaccardIndex(Collection<String> set1,
                                               Collection<String> set2) {
        Set<String> copy1 = new HashSet<String>(set1);
        Set<String> copy2 = new HashSet<String>(set2);
        Set<String> shared = new HashSet<String>(copy1);
        shared.retainAll(copy2);
        copy1.addAll(copy2);
        return (double) shared.size() / copy1.size();
    }

    public static double calculateBinomialPValue(double ratio,
                                                 int sampleSize,
                                                 int success) {
        // Try to use apache math lib
        //BinomialDistribution binomial = distFactory.createBinomialDistribution(sampleSize, ratio);
        //try {
        //return 1.0d - binomial.cumulativeProbability(success - 1);
        //}
        //catch(MathException e) {
        //e.printStackTrace();
        //}
        //return Double.MAX_VALUE;
        // Try to use cern lib
        //One tailed only!!!
        Binomial cernBinomial = new Binomial(sampleSize, ratio, randomeEngine);
        if (success == 0) // To avoid unreasonal value
            success = 1;
        return 1.0d - cernBinomial.cdf(success - 1);
    }

    public static double calOneTailedStandardNormalPvalue(double z) {
        Normal stdNormal = new Normal(0, 1, randomeEngine);
        return stdNormal.cdf(z);
    }

    /**
     * Calculate a up-tailed p-value based on hypergeometic distribution.
     *
     * @param N       the total balls in urn
     * @param s       the white balls (as success)
     * @param n       the sample size
     * @param success the white balls picked up in the sample size
     * @return
     */
    public static double calculateHypergeometricPValue(int N,
                                                       int s,
                                                       int n,
                                                       int success) throws MathException {
        HypergeometricDistributionImpl hyper = new HypergeometricDistributionImpl(N,
                s,
                n);
        return hyper.upperCumulativeProbability(success);
//        int max = Math.min(s, n);
//        if (success >= max / 2)
//            return hyper.upperCumulativeProbability(success);
//        else
//            return 1.0d - hyper.cumulativeProbability(success - 1);
//        return hyper.cumulativeProbability(success, max);
    }

    public static double calTwoTailStandardNormalPvalue(double z) {
        z = Math.abs(z);
        return 2 * (1 - calOneTailedStandardNormalPvalue(z));
    }

    public static double calculateZValue(int success1, int sample1,
                                         int success2, int sample2) {
        double ratio1 = (double) success1 / sample1;
        double ratio2 = (double) success2 / sample2;
        // Need to calculate z value
        double p = (success1 + success2) / (double) (sample1 + sample2);
        double z = (ratio1 - ratio2) / Math.sqrt(p * (1.0 - p) * (1.0 / sample1 + 1.0 / sample2));
        return z;
    }

    public static double twoTailGroupZTest(double p0, double p, int sampleSize) {
        double z = (p - p0) * 1.0 / (Math.sqrt((p0 * (1 - p0)) / sampleSize));
        return calTwoTailStandardNormalPvalue(z);
    }

    public static double chiSquare(double df, double x2) {
        ChiSquare cs = new ChiSquare(df, randomeEngine);
        return cs.cdf(x2);
    }

    public static double calculateEnrichment(double ratio,
                                             int sampleTotal,
                                             int sampleSuccess) {
        double newValue = (double) sampleSuccess / sampleTotal;
        return newValue / ratio;
    }

    public static <T> Set<T> randomSampling(Collection<T> set,
                                            int size) {
        Set<T> rtn = new HashSet<T>();
        int total = set.size();
        List<T> list = null;
        if (set instanceof List)
            list = (List<T>) set;
        else
            list = new ArrayList<T>(set);
        int index;
        while (rtn.size() < size) {
            index = (int) (total * Math.random());
            rtn.add(list.get(index));
        }
        return rtn;
    }

    /**
     * Use this static method to construct a PearsonCorrelation object to get Pearson correaltion
     * value and its p-value.
     *
     * @param values1
     * @param values2
     * @return
     */
    public static PearsonsCorrelation constructPearsonCorrelation(List<Double> values1,
                                                                  List<Double> values2) {
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("Two double lists have different lengths: " +
                    values1.size() + " and " + values2.size());
        Array2DRowRealMatrix matrix = constructMatrix(values1, values2);
        PearsonsCorrelation correlation = new PearsonsCorrelation(matrix);
        return correlation;
    }

    private static Array2DRowRealMatrix constructMatrix(List<Double> values1,
                                                        List<Double> values2) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(values1.size(), 2);
        for (int i = 0; i < values1.size(); i++) {
            matrix.addToEntry(i, 0, values1.get(i));
            matrix.addToEntry(i, 1, values2.get(i));
        }
        return matrix;
    }

    /**
     * Use this method for doing Spearman Correlation.
     *
     * @param values1
     * @param values2
     * @return
     */
    public static SpearmansCorrelation constructSpearmansCorrelation(List<Double> values1,
                                                                     List<Double> values2) {
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("Two double lists have different lengths: " +
                    values1.size() + " and " + values2.size());
        Array2DRowRealMatrix matrix = constructMatrix(values1, values2);
        SpearmansCorrelation correlation = new SpearmansCorrelation(matrix);
        return correlation;
    }

    /**
     * This method is used to calculate Pearson correlation coefficient. The provided two
     * double lists should have the same size.
     *
     * @param values1
     * @param values2
     * @return
     * @deprecate use method constructPearsonCorrealtion instead.
     */
    @Deprecated
    public static Double calculatePearsonCorrelation(List<Double> values1,
                                                     List<Double> values2) {
        // Make sure the two lists have the same size
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("The provided two double lists have different sizes!");
        double av1 = 0.0, av2 = 0.0, var1 = 0.0, var2 = 0.0, var12 = 0.0, c;

        int size = values1.size();
        for (int i = 0; i < size; i++) {
            av1 += values1.get(i);
            av2 += values2.get(i);
        }
        av1 /= (double) size;
        av2 /= (double) size;
        for (int i = 0; i < size; i++) {
            var1 += (values1.get(i) - av1) * (values1.get(i) - av1);
            var2 += (values2.get(i) - av2) * (values2.get(i) - av2);
            var12 += (values1.get(i) - av1) * (values2.get(i) - av2);
        }
        if (var1 * var2 == 0.0) {
            c = 1.0;
        } else {
            c = var12 / Math.sqrt(Math.abs(var1 * var2));
        }

        return c;
    }

    @Deprecated
    public static double calculatePearsonCorrelation1(List<Double> values1,
                                                      List<Double> values2) {
        // Make sure the two lists have the same size
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("The provided two double lists have different sizes!");
        double av1 = 0.0, av2 = 0.0, sqav1 = 0.0, sqav2 = 0.0, cross = 0.0d;

        int size = values1.size();
        for (int i = 0; i < size; i++) {
            av1 += values1.get(i);
            sqav1 += values1.get(i) * values1.get(i);
            av2 += values2.get(i);
            sqav2 += values2.get(i) * values2.get(i);
            cross += values1.get(i) * values2.get(i);
        }
        av1 /= size;
        sqav1 /= size;
        av2 /= size;
        sqav2 /= size;
        cross /= size;
        return (cross - av1 * av2) / Math.sqrt((sqav1 - av1 * av1) * (sqav2 - av2 * av2));
    }

    /**
     * This method should be used instead of randomSampling(Collection<String>, int) since
     * it should be a more standard randomization implementation.
     *
     * @param set
     * @param size
     * @param randomizer
     * @return
     */
    @SuppressWarnings("unchecked")
    public static <T> Set<T> randomSampling(Collection<T> set,
                                            int size,
                                            RandomData randomizer) {
        Object[] objects = randomizer.nextSample(set, size);
        Set<T> rtn = new HashSet<T>();
        for (Object obj : objects)
            rtn.add((T) obj);
        return rtn;
    }

    /**
     * Permutate a list of objects.
     *
     * @param list
     * @param randomizer
     * @return
     */
    public static <T> List<T> permutate(List<T> list,
                                        RandomData randomizer) {
        int[] indices = randomizer.nextPermutation(list.size(), list.size());
        List<T> rtn = new ArrayList<T>();
        for (int i : indices) {
            rtn.add(list.get(i));
        }
        return rtn;
    }

    /**
     * Same as another method permutate() but uses the API from math3.
     *
     * @param list
     * @param randomizer
     * @return
     */
    public static <T> List<T> permutate(List<T> list,
                                        RandomDataGenerator randomizer) {
        int[] indices = randomizer.nextPermutation(list.size(), list.size());
        List<T> rtn = new ArrayList<T>();
        for (int i : indices) {
            rtn.add(list.get(i));
        }
        return rtn;
    }

    /**
     * Random permutate a map
     *
     * @param keyToValue
     * @return
     */
    public static <T, G> Map<T, G> permutate(Map<T, G> keyToValue) {
        Map<T, G> rtn = new HashMap<T, G>();
        List<T> keys = new ArrayList<T>(keyToValue.keySet());
        List<G> values = new ArrayList<G>(keyToValue.values());
        RandomData randomizer = new RandomDataImpl();
        keys = permutate(keys, randomizer);
        values = permutate(values, randomizer);
        for (int i = 0; i < keys.size(); i++) {
            T key = keys.get(i);
            G value = values.get(i);
            rtn.put(key, value);
        }
        return rtn;
    }

    /**
     * Calculate a FDR from a permutation test. The order in an increasing order.
     *
     * @param realValue
     * @param realValues
     * @param randomValues
     * @return
     */
    public static double calculateFDR(double realValue,
                                      List<Double> realValues,
                                      List<Double> randomValues) {
        int realCount = 0;
        for (int i = 0; i < realValues.size(); i++) {
            if (realValues.get(i) > realValue) {
                realCount = i;
                break;
            }
        }
        int randomCount = 0;
        for (int i = 0; i < randomValues.size(); i++) {
            if (randomValues.get(i) > realValue) {
                randomCount = i;
                break;
            }
        }
        if (randomCount == 0)
            return 1.0d;
        double denom = (double) realCount / realValues.size();
        double num = (double) randomCount / randomValues.size();
        double fdr = num / denom;
        if (fdr > 1.0d)
            fdr = 1.0d;
        return fdr;
    }

    /**
     * Use this method to calculate FDR from a list of pvalues using Benjamini-Hochberg
     * method. The implementation of this method is based on the source code for MEMo
     * (http://cbio.mskcc.org/tools/memo/).
     *
     * @param pvalues
     * @return
     */
    public static List<Double> calculateFDRWithBenjaminiHochberg(List<Double> pvalues) {
        // Make sure the parsed list if sorted.
        Collections.sort(pvalues);
        List<Double> fdrs = new ArrayList<Double>(pvalues);
        int size = pvalues.size();
        // The last p-value (biggest) should be the same as FDR.
        for (int i = size - 2; i >= 0; i--) {
            double right = fdrs.get(i + 1);
            double pvalue = pvalues.get(i);
            double left = pvalue * (size / (i + 1));
            double fdr = Math.min(left, right);
            fdrs.set(i, fdr);
        }
        return fdrs;
    }

    public static Double boundDouble01(Double val) {
        return val > 1.0d
                ? 1.0d
                : (val < Double.MIN_NORMAL
                ? Double.MIN_NORMAL
                : val);
    }

    //from https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/blob/master/Python/EmpiricalBrownsMethod.py
    /*
    Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
           A vector of m P-values to combine. May be a list or of type numpy.array.
    Output: A combined P-value using the Empirical Brown's Method (EBM).
            If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
     */
    public static Double calculatePValuesWithEmpiricalBrownsMethod(List<List<Double>> matrix, List<Double> pvalues) {
        RealMatrix covarianceMatrix = calculateCovariances(matrix);
        return combinePValuesWithEBM(covarianceMatrix, pvalues);
    }

    private static Double combinePValuesWithEBM(RealMatrix covarianceMatrix, List<Double> pvalues) {
        Integer m = covarianceMatrix.getRowDimension();
        Double df_fisher = 2.0d * m;
        Double Expected = 2.0d * m;
        Double cov_sum = 0.0d;
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < m; j++) {
                cov_sum += covarianceMatrix.getEntry(i, j);
            }
        }
        Double Var = 4.0d * m + 2.0d * cov_sum;
        Double c = Var / (2.0d * Expected);
        Double df_brown = 2.0d * Math.pow(Expected, 2.0d) / Var;
        if (df_brown > df_fisher) {
            df_brown = df_fisher;
            c = 1.0d;
        }
        Double sum_nlog_pvalues = 0.0;
        for (Double pvalue : pvalues) {
            sum_nlog_pvalues -= Math.log(pvalue);
        }
        Double x = 2.0d * sum_nlog_pvalues;

        ChiSquaredDistribution chiSquaredDistributionBrown = new ChiSquaredDistribution(df_brown);
        return 1.0d - chiSquaredDistributionBrown.cumulativeProbability(1.0d * x / c);
    }

    private static double[][] transformMatrix(List<List<Double>> matrix) {
        double[][] transformedMatrix = new double[matrix.size()][matrix.get(0).size()];
        for (int i = 0; i < matrix.size(); i++) {
            List<Double> vector = matrix.get(i);
            double[] transformedVector = new double[vector.size()];
            double[] sVector = new double[vector.size()];
            SummaryStatistics vectorSummary = new SummaryStatistics();
            for (Double element : vector) {
                vectorSummary.addValue(element);
            }
            Double mean = vectorSummary.getMean();
            Double std = vectorSummary.getStandardDeviation();
            for (int j = 0; j < vector.size(); j++) {
                Double element = vector.get(j);
                sVector[j] = (element - mean) / std;
            }
            for (int k = 0; k < sVector.length; k++) {
                transformedVector[k] = -2.0d * Math.log(calculateEmpiricalCDF(sVector, sVector[k]));
            }
            transformedMatrix[i] = transformedVector;
        }
        return transformedMatrix;
    }

    private static Double calculateEmpiricalCDF(double[] vector, Double value) {
        Integer nltev = 0;
        for (double element : vector) {
            nltev += element <= value
                    ? 1
                    : 0;
        }
        return (double) nltev / (double) vector.length;
    }

    private static RealMatrix calculateCovariances(List<List<Double>> matrix) {
        double[][] transformedMatrix = transformMatrix(matrix);
        double[][] columnwiseTransformedMatrix = convertToColumnwise(transformedMatrix);
        Covariance covarianceMatrix = new Covariance(columnwiseTransformedMatrix);
        return covarianceMatrix.getCovarianceMatrix();
    }

    private static double[][] convertToColumnwise(double[][] rowwise) {
        double[][] columnwise = new double[rowwise[0].length][rowwise.length];
        for (int rowIdx = 0; rowIdx < rowwise.length; rowIdx++) {
            for (int colIdx = 0; colIdx < rowwise[0].length; colIdx++) {
                columnwise[colIdx][rowIdx] = rowwise[rowIdx][colIdx];
            }
        }
        return columnwise;
    }

    /**
     * This method is used to combine a collection of pvalues using Fisher's method.
     * (see http://en.wikipedia.org/wiki/Fisher's_method).
     *
     * @param pvalues
     */
    public static double combinePValuesWithFisherMethod(Collection<Double> pvalues) throws MathException {
        // Have to make sure there is no zero in the pvalues collection. Otherwise,
        // log will throw an exception
        for (Double pvalue : pvalues) {
            if (pvalue.equals(0.0d)) {
                throw new IllegalArgumentException("The pvalue list contains 0, which is not allowed!");
            }
        }
        Double total = 0.0d;
        for (Double pvalue : pvalues) {
            total += Math.log(pvalue);
        }
        total *= -2.0d;
        ChiSquaredDistribution distribution = new ChiSquaredDistribution(2 * pvalues.size());
        return 1.0d - distribution.cumulativeProbability(total);
    }

    /**
     * This implementation of proportion test is simplified based on the R prop.test. The chisq correction
     * is based on book: Introductory to Statistics with R.
     */
    public static double proportionTest(int success,
                                        int sampleSize,
                                        double proportation) throws MathException {
        double diff = Math.abs(success - sampleSize * proportation);
        // Check what correction should be used
        double yatesCorrection = Math.min(0.50d, diff);
        double chisqr = (diff - yatesCorrection) * (diff - yatesCorrection) / (sampleSize * proportation * (1.0d - proportation));
        // Calculate chisq pvalues
        ChiSquaredDistribution distribution = new ChiSquaredDistribution(1);
        return 1.0d - distribution.cumulativeProbability(chisqr);
    }

    @Test
    public void test() throws Exception {
//        double p = calculateHypergeometricPValue(45,
//                                                 17,
//                                                 15,
//                                                 9);
//        System.out.println("P value: " + p);
        double zvalue = calculateZValue(76, 91, 14, 18);
        double pvalue = calTwoTailStandardNormalPvalue(zvalue);
        System.out.println(pvalue);
    }

    @Test
    public void testHyper() throws MathException {
        HypergeometricDistribution hyper = new HypergeometricDistributionImpl(9083, 1977, 201);
        double pvalue1 = hyper.cumulativeProbability(54, 201);
        System.out.println(pvalue1);
        double pvalue2 = hyper.cumulativeProbability(53);
        System.out.println(pvalue2);
        System.out.println(pvalue1 + pvalue2);
    }

    public Normal normalPDF() {
        return new Normal(0, 1, randomeEngine);
    }

    @Test
    public void testProportionTest() throws MathException {
        Double pvalue = proportionTest(16, 560, 0.012);
        System.out.println("Pvalue for 16, 560, 0.012: " + pvalue);
        pvalue = proportionTest(84, 560, 0.012);
        System.out.println("Pvalue for 84, 560, 0.012: " + pvalue);
        pvalue = proportionTest(6, 560, 0.012);
        System.out.println("Pvalue for 6, 560, 0.012: " + pvalue);
        pvalue = proportionTest(13, 560, 0.016);
        System.out.println("Pvalue for 13, 560, 0.016: " + pvalue);
    }

    @Test
    public void testRandomization() {
        int total = 500;
        RandomDataImpl randomizer = new RandomDataImpl();
        for (int i = 0; i < 10; i++) {
            StringBuilder builder = new StringBuilder();
            for (int j = 0; j < 100; j++) {
                //int random = (int) (total * Math.random());
                int random = randomizer.nextInt(0, total - 1);
                builder.append(random).append(", ");
            }
            System.out.println(builder.toString());
        }
    }

    @Test
    public void testTestdata() {
        int rowLen = row1.length;
        Assert.assertEquals(rowLen, row2.length);
        Assert.assertEquals(rowLen, row3.length);
        Assert.assertEquals(rowLen, row4.length);
        Assert.assertEquals(rowLen, row5.length);
        Assert.assertEquals(rowLen, row6.length);
        Assert.assertEquals(rowLen, row7.length);
    }

    @Test
    public void testCalculatePValuesWithEmpiricalBrownsMethod() {
        List<Double> row1 = new ArrayList<>(Arrays.asList(this.row1));
        List<Double> row2 = new ArrayList<>(Arrays.asList(this.row2));
        List<Double> row3 = new ArrayList<>(Arrays.asList(this.row3));
        List<Double> row4 = new ArrayList<>(Arrays.asList(this.row4));
        List<Double> row5 = new ArrayList<>(Arrays.asList(this.row5));
        List<Double> row6 = new ArrayList<>(Arrays.asList(this.row6));
        List<Double> row7 = new ArrayList<>(Arrays.asList(this.row7));

        List<List<Double>> matrix = new ArrayList<>();
        matrix.add(row1);
        matrix.add(row2);
        matrix.add(row3);
        matrix.add(row4);
        matrix.add(row5);
        matrix.add(row6);
        matrix.add(row7);

        List<Double> pvalues = new ArrayList<>(Arrays.asList(this.glypPvals));

        Double actualBrowns = calculatePValuesWithEmpiricalBrownsMethod(matrix, pvalues);

        Double expectedBrowns = 4.821679e-07;
        Assert.assertEquals("Browns method p-value calculation incorrect", expectedBrowns, actualBrowns, 1.0e-13);
    }

    @Test
    public void testCombinePValuesWithFisherMethod() throws MathException {
        List<Double> pvalues = new ArrayList<>(Arrays.asList(this.glypPvals));

        Double actualFishers = combinePValuesWithFisherMethod(pvalues);

        Double expectedFishers = 1.438732e-08;
        Assert.assertEquals("Fishers method p-value calculation incorrect", expectedFishers, actualFishers, 1.0e-13);
    }
}
