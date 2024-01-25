package OSM;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import OSM.XMLParser.Element;
import OSM.XMLParser.InnerNodeListener;
import OSM.XMLParser.XMLNodeListener;

import static javax.xml.stream.XMLStreamConstants.*;
import static OSM.Util.*;

public abstract class OSMReader {

    protected File osm;

    public OSMReader(File osm) {
        this.osm = osm;
    }
    public abstract List<Way> readWays();
    public abstract List<Node> readNodes(final Map<Long, Node> streetNodes);

    public static class Way {
        public Way(List<Long> nodeIDs, String streetType) {
            this.nodeIDs = nodeIDs;
            this.streetType = streetType;
            this.network = "";
            this.operator = "";
        }
        public Way(List<Long> nodeIDs, String network, String operator) {
            this.nodeIDs = nodeIDs;
            this.streetType = "";
            this.network = network;
            this.operator = operator;
        }
        public List<Long> nodeIDs;
        public String streetType;
        public String network;
        public String operator;
    }

    public static class StreamReader extends OSMReader {
        protected XMLInputFactory factory = XMLInputFactory.newInstance();
        public StreamReader(File osm) {
            super(osm);
        }
        @Override
        public List<Way> readWays() {
            List<Way> ways = new ArrayList<Way>();
            try {
                FileReader fr = new FileReader(osm);
                XMLStreamReader wayReader = factory.createXMLStreamReader(fr);
                while (wayReader.hasNext()) {
                    wayReader.next();
                    if (wayReader.getEventType() == START_ELEMENT && wayReader.getName().toString().equals("way")) {
                        String area = "no";
                        String highway = "";
                        String tracktype = "";
                        boolean closed = false;
                        List<Long> nd = new ArrayList<Long>();
                        Long id = new Long(getAttribute(wayReader, "id", "-1"));
                        while (wayReader.hasNext() && !closed) {
                            wayReader.next();
                            if (wayReader.getEventType() == START_ELEMENT && wayReader.getName().toString().equals("tag")) {
                                if (getAttribute(wayReader, "k").equals("area")) {
                                    area = getAttribute(wayReader, "v", area);
                                }
                                if (getAttribute(wayReader, "k").equals("highway")) {
                                    highway = getAttribute(wayReader, "v", highway);
                                }
                                if (getAttribute(wayReader, "k").equals("tracktype")) {
                                    tracktype = getAttribute(wayReader, "v", tracktype);
                                }
                            }
                            if (wayReader.getEventType() == START_ELEMENT && wayReader.getName().toString().equals("nd")) {
                                nd.add(new Long(getAttribute(wayReader, "ref")));
                            }
                            if (wayReader.getEventType() == END_ELEMENT && wayReader.getName().toString().equals("way")) {
                                closed = true;
                            }
                        }
                        if (streetTypes.containsKey(highway) && !area.equals("yes")) {
                            if (streetTypes.containsKey(highway + ":" + tracktype)) {
                                highway = highway + ":" + tracktype;
                            }
                            ways.add(new Way(nd, highway));
                        }
                    }
                }
                wayReader.close();
                fr.close();
            } catch (NumberFormatException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (XMLStreamException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return ways;
        }
        @Override
        public List<Node> readNodes(Map<Long, Node> streetNodes) {
            List<Node> nodes = new ArrayList<Node>();
            try {
                FileReader fr = new FileReader(osm);
                XMLStreamReader nodeReader = factory.createXMLStreamReader(fr);
                while (nodeReader.hasNext()) {
                    nodeReader.next();
                    if (nodeReader.getEventType() == START_ELEMENT && nodeReader.getName().toString().equals("node")) {
                        Long id = new Long(getAttribute(nodeReader, "id", "-1"));
                        double lat = Double.parseDouble(getAttribute(nodeReader, "lat", "0"));
                        double lon = Double.parseDouble(getAttribute(nodeReader, "lon", "0"));
                        if (streetNodes.containsKey(id)) {
                            Node node = new Node(id, lat, lon);
                            streetNodes.put(id, node);
                            nodes.add(node);
                        }
                    }
                }
                nodeReader.close();
                fr.close();
            } catch (NumberFormatException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (XMLStreamException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return nodes;
        }
        public static String getAttribute(XMLStreamReader reader, String key) {
            return getAttribute(reader, key, "");
        }
        public static String getAttribute(XMLStreamReader reader, String key, String defaultValue) {
            for (int i = 0; i < reader.getAttributeCount(); i++) {
                if (reader.getAttributeName(i).toString().equals(key)) {
                    defaultValue = reader.getAttributeValue(i);
                }
            }
            return defaultValue;
        }
    }

    public static class StreamBikeSharingReader extends OSMReader {
        protected XMLInputFactory factory = XMLInputFactory.newInstance();
        public StreamBikeSharingReader(File osm) {
            super(osm);
        }
        @Override
        public List<Way> readWays() {
            List<Way> ways = new ArrayList<Way>();
            try {
                FileReader fr = new FileReader(osm);
                XMLStreamReader wayReader = factory.createXMLStreamReader(fr);
                while (wayReader.hasNext()) {
                    wayReader.next();
                    if (wayReader.getEventType() == START_ELEMENT && wayReader.getName().toString().equals("way")) {
                        String network = "";
                        String operator = "";
                        boolean bicycleRntal = false;
                        boolean closed = false;
                        List<Long> nd = new ArrayList<Long>();
                        Long id = new Long(getAttribute(wayReader, "id", "-1"));
                        while (wayReader.hasNext() && !closed) {
                            wayReader.next();
                            if (wayReader.getEventType() == START_ELEMENT && wayReader.getName().toString().equals("tag")) {
                                if (getAttribute(wayReader, "k").equals("amenity") && getAttribute(wayReader, "v").equals("bicycle_rental")) {
                                    bicycleRntal = true;
                                }
                                if (getAttribute(wayReader, "k").equals("network")) {
                                    network = getAttribute(wayReader, "v", network);
                                }
                                if (getAttribute(wayReader, "k").equals("operator")) {
                                    operator = getAttribute(wayReader, "v", operator);
                                }
                            }
                            if (wayReader.getEventType() == START_ELEMENT && wayReader.getName().toString().equals("nd")) {
                                nd.add(new Long(getAttribute(wayReader, "ref")));
                            }
                            if (wayReader.getEventType() == END_ELEMENT && wayReader.getName().toString().equals("way")) {
                                closed = true;
                            }
                        }
                        if (bicycleRntal) {
                            ways.add(new Way(nd, network, operator));
                        }
                    }
                }
                wayReader.close();
                fr.close();
            } catch (NumberFormatException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (XMLStreamException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return ways;
        }
        @Override
        public List<Node> readNodes(Map<Long, Node> wayNodes) {
            List<Node> nodes = new ArrayList<Node>();
            try {
                FileReader fr = new FileReader(osm);
                XMLStreamReader nodeReader = factory.createXMLStreamReader(fr);
                while (nodeReader.hasNext()) {
                    nodeReader.next();
                    if (nodeReader.getEventType() == START_ELEMENT && nodeReader.getName().toString().equals("node")) {
                        String network = "";
                        String operator = "";
                        boolean bicycleRntal = false;
                        boolean closed = false;
                        Long id = new Long(getAttribute(nodeReader, "id", "-1"));
                        double lat = Double.parseDouble(getAttribute(nodeReader, "lat", "0"));
                        double lon = Double.parseDouble(getAttribute(nodeReader, "lon", "0"));
                        while (nodeReader.hasNext() && !closed) {
                            nodeReader.next();
                            if (nodeReader.getEventType() == START_ELEMENT && nodeReader.getName().toString().equals("tag")) {
                                if (getAttribute(nodeReader, "k").equals("amenity") && getAttribute(nodeReader, "v").equals("bicycle_rental")) {
                                    bicycleRntal = true;
                                }
                                if (getAttribute(nodeReader, "k").equals("network")) {
                                    network = getAttribute(nodeReader, "v", network);
                                }
                                if (getAttribute(nodeReader, "k").equals("operator")) {
                                    operator = getAttribute(nodeReader, "v", operator);
                                }
                            }
                            if (nodeReader.getEventType() == END_ELEMENT && nodeReader.getName().toString().equals("node")) {
                                closed = true;
                            }
                        }
                        if (bicycleRntal) {
                            NodeBikeSharing node = new NodeBikeSharing(id, lat, lon, network, operator);
                            nodes.add(node);
                        }
                        if (wayNodes.containsKey(id)) {
                            Node node = new Node(id, lat, lon);
                            wayNodes.put(id, node);
                        }
                    }
                }
                nodeReader.close();
                fr.close();
            } catch (NumberFormatException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (XMLStreamException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return nodes;
        }
        public static String getAttribute(XMLStreamReader reader, String key) {
            return getAttribute(reader, key, "");
        }
        public static String getAttribute(XMLStreamReader reader, String key, String defaultValue) {
            for (int i = 0; i < reader.getAttributeCount(); i++) {
                if (reader.getAttributeName(i).toString().equals(key)) {
                    defaultValue = reader.getAttributeValue(i);
                }
            }
            return defaultValue;
        }
    }

    public static class XMLParserReader extends OSMReader {
        public XMLParserReader(File osm) {
            super(osm);
        }
        @Override
        public List<Way> readWays() {
            final List<Way> ways = new ArrayList<Way>();
            XMLParser parser = new XMLParser();
            parser.addNodeListener(new XMLNodeListener() {
                private String area = "no";
                private String highway = "";
                private String tracktype = "";
                private List<Long> nd = new ArrayList<Long>();
                private Long id = null;
                private Map<String, XMLNodeListener> subNodeListeners = new HashMap<String, XMLNodeListener>(); {
                    subNodeListeners.put("tag", new InnerNodeListener("tag") {
                        public void nodeBegin(Element e) {
                            if (e.hasAttribute("k", "area")) {
                                area = e.getAttribute("v", area);
                            }
                            if (e.hasAttribute("k", "highway")) {
                                highway = e.getAttribute("v", highway);
                            }
                            if (e.hasAttribute("k", "tracktype")) {
                                tracktype = e.getAttribute("v", tracktype);
                            }
                        }
                    });
                    subNodeListeners.put("nd", new InnerNodeListener("nd") {
                        public void nodeBegin(Element e) {
                            nd.add(new Long(e.getAttribute("ref")));
                        }
                    });
                }
                public String getNodeName() {return "way";}
                public Map<String, XMLNodeListener> getSubNodeListeners() {return subNodeListeners;}
                public void nodeBegin(Element e) {
                    id = new Long(e.getAttribute("id", "-1"));
                }
                public void nodeEnd() {
                    if (streetTypes.containsKey(highway) && !area.equals("yes")) {
                        if (streetTypes.containsKey(highway + ":" + tracktype)) {
                            highway = highway + ":" + tracktype;
                        }
                        ways.add(new Way(nd, highway));
                    }
                }
            });
            parser.process(osm);
            return ways;
        }
        @Override
        public List<Node> readNodes(final Map<Long, Node> streetNodes) {
            final List<Node> nodes = new ArrayList<Node>();
            XMLParser parser = new XMLParser();
            parser.addNodeListener(new InnerNodeListener("node") {
                public void nodeBegin(Element e) {
                    Long id = new Long(e.getAttribute("id", "-1"));
                    double lat = Double.parseDouble(e.getAttribute("lat", "0"));
                    double lon = Double.parseDouble(e.getAttribute("lon", "0"));
                    if (streetNodes.containsKey(id)) {
                        Node node = new Node(id, lat, lon);
                        streetNodes.put(id, node);
                        nodes.add(node);
                    }
                }
            });
            parser.process(osm);
            return nodes;
        }
    }

}
