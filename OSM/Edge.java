package OSM;

import static OSM.Util.*;

public class Edge implements Comparable<Edge> {

    private Node from;
    private Node to;
    private String streetType;

    public Edge(Node from, Node to, String streetType) {
        this.from = from;
        this.to = to;
        this.streetType = streetType;
    }

    @Override
    public int compareTo(Edge o) {
        if (this.getMinNode().getOSMId() == o.getMinNode().getOSMId()) {
            return (int) Math.signum(this.getMaxNode().getOSMId() - o.getMaxNode().getOSMId());
        } else {
            return (int) Math.signum(this.getMinNode().getOSMId() - o.getMinNode().getOSMId());
        }
    }
    
    @Override
    public boolean equals(Object o) {
        if (o != null && o instanceof Edge) {
            Edge e = (Edge) o;
            return (this.getMaxNode() == e.getMaxNode()) && (this.getMinNode() == e.getMinNode());
        }
        return false;
    }

    @Override
    public String toString() {
        return "a " + (from.getGraphId() + 1) + " " + (to.getGraphId() + 1) + " " + getTravelTime();
    }
    
    public double getDistance() {
        return this.from.getDist(this.to);
    }
    
    public double getSpeed() {
        return streetTypes.get(this.getStreetType());
    }
    
    public int getTravelTime() {
        return (int) Math.max(1.0, Math.ceil(((getDistance() / getSpeed()) * 0.036)));
    }

    public Node getMinNode() {
        return this.getFrom().getOSMId() < this.getTo().getOSMId() ? this.getFrom() : this.getTo();
    }

    public Node getMaxNode() {
        return this.getFrom().getOSMId() < this.getTo().getOSMId() ? this.getTo() : this.getFrom();
    }

    public Node getFrom() {
        return from;
    }

    public void setFrom(Node from) {
        this.from = from;
    }

    public Node getTo() {
        return to;
    }

    public void setTo(Node to) {
        this.to = to;
    }

    public String getStreetType() {
        return streetType;
    }

    public void setStreetType(String streetType) {
        this.streetType = streetType;
    }
    
}
