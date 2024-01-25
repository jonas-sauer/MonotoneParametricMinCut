package OSM;

import static OSM.Util.*;

public class NodeBikeSharing extends Node {

    private String network;

    private String operator;

    public NodeBikeSharing(long OSMId, double lat, double lon, String network, String operator) {
        super(OSMId, lat, lon);
        this.network = network;
        this.operator = operator;
    }

    @Override
    public String toString() {
        return super.toString() + " " + this.getNetwork() + " " + this.getOperator();
    }

    public String getNetwork() {
        return this.network;
    }

    public void setNetwork(String network) {
        this.network = network;
    }

    public String getOperator() {
        return this.operator;
    }

    public void setOperator(String operator) {
        this.operator = operator;
    }

}
