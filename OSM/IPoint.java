package OSM;

import java.util.Comparator;

public interface IPoint<T> {
    
    public double getX();    
    public double getY();
    
    public double getDistX(T p);    
    public double getDistY(T p);    
    public double getDist(T p);
    
    public static class CompareByY<T extends IPoint<? super T>> implements Comparator<T> {
        public int compare(T a, T b) {
            return (int) Math.signum(a.getY() - b.getY());
        }
    }
    
    public static class CompareByX<T extends IPoint<? super T>> implements Comparator<T> {
        public int compare(T a, T b) {
            return (int) Math.signum(a.getX() - b.getX());
        }
    }

}
