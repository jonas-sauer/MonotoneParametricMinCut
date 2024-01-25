package OSM;

import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import static OSM.Util.*;
import static OSM.OSMReader.*;

public class Parser {

    public static class Graph {
        public Graph(Node[] nodes, List<Edge> edges) {
            this.nodes = nodes;
            this.edges = edges;
        }
        public Node[] nodes;
        public List<Edge> edges;
    }

    private Parser() {}

    public static Graph parseOSM(File osm) {
        long time = System.nanoTime();
        OSMReader reader = new StreamReader(osm);

        // read ways
        System.out.println("Reading ways");
        List<Way> ways = reader.readWays();
        Map<Long, Node> streetNodes = new HashMap<Long, Node>();
        for (Way way : ways) {
            for (Long nodeID : way.nodeIDs) {
                streetNodes.put(nodeID, null);
            }
        }
        System.out.println("\tmaxNumNodes: " + streetNodes.size());

        // read nodes
        System.out.println("Reading nodes");
        List<Node> vertices = reader.readNodes(streetNodes);
        System.out.println("\tmaxNumNodes: " + vertices.size());

        // contract nearby nodes
        System.out.println("Contracting nearby nodes");
        PointReduction<Node> pointReduction = new PointReduction<Node>(vertices, MIN_NODE_DIST_IN_CENTIMETRE);
        System.out.println("\tnumNodeReductions: " + pointReduction.getNumReductions());

        // remove incomplete edges & nodes
        System.out.println("Removing incomplete edges and nodes");
        vertices.clear();
        List<Edge> edges = new LinkedList<Edge>();
        for (Way way : ways) {
            for (int i = 1; i < way.nodeIDs.size(); i++) {
                Node from = streetNodes.get(way.nodeIDs.get(i - 1));
                Node to = streetNodes.get(way.nodeIDs.get(i));
                if (from != null && to != null) {
                    from = pointReduction.getParent(from);
                    to = pointReduction.getParent(to);
                    edges.add(new Edge(from, to, way.streetType));
                    if (!from.isUsed()) {
                        vertices.add(from);
                        from.setUsed(true);
                    }
                    if (!to.isUsed()) {
                        vertices.add(to);
                        to.setUsed(true);
                    }
                }
            }
        }
        System.out.println("\tmaxNumEdges: " + edges.size());
        System.out.println("\tnumNodes: " + vertices.size());

        // remove double edges & loops
        System.out.println("Removing double edges and loops");
        Collections.sort(edges);
        Edge last = null;
        for (ListIterator<Edge> i = edges.listIterator(); i.hasNext();) {
            Edge next = i.next();
            if (next.getFrom() == next.getTo()) {
                i.remove();
            } else if (next.equals(last)) {
                i.remove();
            } else {
                last = next;
            }
        }
        System.out.println("\tnumEdges: " + edges.size());

        // Create Output
        System.out.println("Creating output");
        Node[] nodes = new Node[vertices.size()];
        int i = 0;
        for (Node n : vertices) {
            n.setGraphId(i);
            nodes[i++] = n;
        }
        time = System.nanoTime() - time;
        System.out.println("\telapsedTime: " + ((time / 1000) / 1000.0) + "ms");

        return new Graph(nodes, edges);
    }

    public static List<Node> parseBickeSharing(File osm) {
        long time = System.nanoTime();
        OSMReader reader = new StreamBikeSharingReader(osm);

        // read ways
        System.out.println("Reading ways");
        List<Way> ways = reader.readWays();
        Map<Long, Node> streetNodes = new HashMap<Long, Node>();
        for (Way way : ways) {
            for (Long nodeID : way.nodeIDs) {
                streetNodes.put(nodeID, null);
            }
        }
        System.out.println("\tnumNodes: " + ways.size());

        // read nodes
        System.out.println("Reading nodes");
        List<Node> vertices = reader.readNodes(streetNodes);
        System.out.println("\tmaxNumNodes: " + (vertices.size() + ways.size()));

        // remove incomplete edges & nodes
        System.out.println("Creating Nodes from ways");
        List<Edge> edges = new LinkedList<Edge>();
        for (Way way : ways) {
            GeographicPoint min = new GeographicPoint(360, 360);
            GeographicPoint max = new GeographicPoint(-360, -360);
            for (int i = 0; i < way.nodeIDs.size(); i++) {
                Node node = streetNodes.get(way.nodeIDs.get(i));
                if (node != null) {
                    min = min.minimize(node);
                    max = max.maximize(node);
                }
            }
            if (min.getLat() < 360) {
                NodeBikeSharing node = new NodeBikeSharing(0, (min.getLat() + max.getLat()) / 2.0, (min.getLon() + max.getLon()) / 2.0, way.network, way.operator);
                vertices.add(node);
            }
        }
        System.out.println("\tnumNodes: " + vertices.size());

        return vertices;
    }

}
