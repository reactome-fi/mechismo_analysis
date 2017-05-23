/*
 * Created on May 22, 2017
 *
 */
package org.reactome.r3.util;

import java.awt.Dimension;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.ui.ApplicationFrame;

/**
 * Plotter based on JFreeChart API
 * @author gwu
 *
 */
public class Plotter extends ApplicationFrame {
    
//    static {
//       ChartFactory.setChartTheme(StandardChartTheme.createJFreeTheme());
//    }
    
    /**
     * Default constructor.
     */
    public Plotter(String title) {
        super(title);
    }
    
    /**
     * Plot multiple datasets as a histogram.
     * @param dataSetToValues
     * @param bins
     * @param mainTitle
     * @param xLabel
     * @param yLabel
     */
    public void plotHistograpm(Map<String, List<Double>> dataSetToValues,
                               int bins,
                               String mainTitle,
                               String xLabel,
                               String yLabel) {
        // Get the maximum and minimum values
        double[] minMax = {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        dataSetToValues.values()
        .stream()
        .flatMap(List::stream)
        .distinct()
        .forEach(value -> {
            minMax[0] = Math.min(minMax[0], value);
            minMax[1] = Math.max(minMax[1], value);
        });
        HistogramDataset dataset = new HistogramDataset();
        dataSetToValues.keySet().stream()
                                .sorted((set1, set2) -> {
                                    List<Double> values1 = dataSetToValues.get(set1);
                                    List<Double> values2 = dataSetToValues.get(set2);
                                    return values1.size() - values2.size(); 
                                 })
                                .forEach(name -> {
                                    List<Double> values = dataSetToValues.get(name);
                                    double[] aValues = values.stream().mapToDouble(v -> v.doubleValue()).toArray();
                                    dataset.addSeries(name,
                                                      aValues,
                                                      bins,
                                                      minMax[0],
                                                      minMax[1]);
                                });
        // create the chart...
        JFreeChart chart1 = ChartFactory.createHistogram(mainTitle,
                                                        xLabel, 
                                                        yLabel, 
                                                        dataset,
                                                        PlotOrientation.VERTICAL, 
                                                        true, 
                                                        true, 
                                                        true);
        JFreeChart chart = chart1;
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setFillZoomRectangle(true);
//        chartPanel.setMouseWheelEnabled(true);
        chartPanel.setPreferredSize(new Dimension(1000, 540));
        setContentPane(chartPanel);
        pack();
        setVisible(true);
    }
    
}
