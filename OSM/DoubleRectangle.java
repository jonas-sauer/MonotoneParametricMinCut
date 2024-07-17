package OSM;

public class DoubleRectangle {

    private double minX;
    private double minY;
    private double maxX;
    private double maxY;

    public DoubleRectangle(double minX, double minY, double maxX, double maxY) {
        this.minX = minX;
        this.minY = minY;
        this.maxX = maxX;
        this.maxY = maxY;
    }

    public boolean contains(double x, double y) {
        return (minX <= x) && (x <= maxX) && (minY <= y) && (y <= maxY); 
    }

    public double getMinX() {
        return minX;
    }

    public double getMinY() {
        return minY;
    }

    public double getMaxX() {
        return maxX;
    }

    public double getMaxY() {
        return maxY;
    }

}
