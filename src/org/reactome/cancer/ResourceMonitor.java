package org.reactome.cancer;

public class ResourceMonitor {
    private Long startMethodTime;
    private Long startLoopTime;
    private Double maxMemUsed;

    public ResourceMonitor() {
        this.maxMemUsed = CalculateJavaMemFootprintGiB();
    }

    public void StartMethodTimer() {
        this.startMethodTime = System.currentTimeMillis();
    }

    public void EndMethodTimer() {
        Long endMethodTime = System.currentTimeMillis();
        System.out.println(String.format("Completed after %.2f minutes\n" +
                        "Max memory footprint: %.2f GiB",
                (endMethodTime - this.startMethodTime) / 60000.0,
                this.maxMemUsed));
    }

    public void StartLoopTimer() {
        this.startLoopTime = System.currentTimeMillis();
    }

    public void EndLoopTimer(String iterationName) {
        Long endLoopTime = System.currentTimeMillis();
        System.out.println(String.format("Completed permutation iteration %s in %.2f minutes\n" +
                        "Total running time: %.2f minutes\n" +
                        "Max memory footprint: %.2f GiB",
                iterationName,
                (endLoopTime - this.startLoopTime) / 60000.0,
                (endLoopTime - this.startMethodTime) / 60000.0,
                this.maxMemUsed));
    }

    public void CalculateMemUsed() {
        double curMemUsed = CalculateJavaMemFootprintGiB();
        this.maxMemUsed = curMemUsed > this.maxMemUsed ?
                curMemUsed :
                this.maxMemUsed;
    }

    private double CalculateJavaMemFootprint() {
        return Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
    }

    private double CalculateJavaMemFootprintGiB() {
        return CalculateJavaMemFootprint() / Math.pow(2.0, 30);
    }
}
