package OSM;

import static OSM.Util.*;

public class Node extends GeographicPoint {

    private final long OSMId;
        
    private int GraphId;
    
    private boolean used;

    public Node(long OSMId, double lat, double lon) {
        super(lat, lon);
        this.OSMId = OSMId;
        this.setUsed(false);
    }

    @Override
    public String toString() {
        return "v " + (this.getGraphId() + 1) + " " + ((int) (this.getLon() * 1000000)) + " " + ((int) (this.getLat() * 1000000));
    }

    public int getGraphId() {
        return this.GraphId;
    }

    public void setGraphId(int graphId) {
        this.GraphId = graphId;
    }

    public long getOSMId() {
        return this.OSMId;
    }

    public boolean isUsed() {
        return this.used;
    }

    public void setUsed(boolean used) {
        this.used = used;
    }

}
