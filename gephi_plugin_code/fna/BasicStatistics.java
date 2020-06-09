package org.msk.supernodehierarchy;

public class BasicStatistics {
    Double[] data;
    int size;   

    public BasicStatistics(Double[] data) {
        this.data = data;
        size = data.length;
    }   

    Double getMean() {
        Double sum = 0.0;
        for(Double a : data)
            sum += a;
        return sum/size;
    }

    Double getVariance() {
        Double mean = getMean();
        Double temp = 0.0;
        for(Double a :data)
            temp += (a-mean)*(a-mean);
        return temp/(size-1);
    }

    Double getStdDev() {
        return Math.sqrt(getVariance());
    }
}


