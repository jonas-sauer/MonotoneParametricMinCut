package OSM;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static OSM.IPoint.*;

public class PointReduction<T extends IPoint<? super T>> {
    
    public static final double COLUMN_WIDTH = 10;

    protected List<T> points; 
    protected Map<T, Label> labels;
    
    protected double minDist;
    protected int numReductions = 0;
    
    public PointReduction(List<T> points, double minDist) {
        this.points = points;
        this.minDist = minDist;
        this.makeLabels();
        this.reduce();
    }
    
    public T getParent(T point) {
        return this.labels.get(point).parent;
    }
    
    protected void makeLabels() {
        this.labels = new HashMap<T, Label>();
        for (T p : this.points) {
            this.labels.put(p, new Label(p));
        }
    }

    protected void reduce() {
        List<List<T>> columns = makeColumns();
        int[] columnIndices = new int[columns.size()];
        Collections.sort(this.points, new CompareByX<T>());
        for (T p : this.points) {
            Label label = this.labels.get(p);
            columnIndices[label.column] = columnIndices[label.column] + 1;
            if (label.parent == p) {
                for (int offset = -1; offset <= 1; offset++) {
                    boolean done = false;
                    List<T> column = columns.get(label.column + offset);
                    for (int j = columnIndices[label.column + offset]; j < column.size() && !done; j++) {
                        T test = column.get(j);
                        if (p.getDist(test) < this.minDist) {
                            this.labels.get(test).parent = p;
                            this.numReductions++;
                        } else {
                            done = p.getDistX(test) > this.minDist;
                        }
                    }
                }
            }
        }
    }
    
    protected List<List<T>> makeColumns () {
        Collections.sort(this.points, new CompareByY<T>());
        List<List<T>> columns = new ArrayList<List<T>>();
        columns.add(new ArrayList<T>());
        List<T> column = new ArrayList<T>();
        T first = this.points.get(0);
        for (T p : this.points) {
            if (first.getDistY(p) > COLUMN_WIDTH * this.minDist) {
                Collections.sort(column, new CompareByX<T>());
                columns.add(column);
                column = new ArrayList<T>();
                first = p;                
            }
            column.add(p);
            this.labels.get(p).column = columns.size();
        }
        Collections.sort(column, new CompareByX<T>());
        columns.add(column);
        columns.add(new ArrayList<T>());
        return columns;
    }
    
    private class Label {
        public T parent;
        public int column;
        public Label(T parent) {
            this.parent = parent;
            this.column = -1;
        }
    }

    public int getNumReductions() {
        return numReductions;
    }
    
}
