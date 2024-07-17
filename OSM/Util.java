package OSM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

public class Util {

    public static final int MIN_NODE_DIST_IN_CENTIMETRE = 100;

    public static final int OPEN = 0;
    public static final int SAVE = 1;

    public static HashMap<String, Double> streetTypes;
    static {
    	streetTypes = new HashMap<String, Double>();
    	streetTypes.put("crossing",                4.5);
    	streetTypes.put("footway",                 4.5);
    	streetTypes.put("living_street",           4.5);
    	streetTypes.put("motorway",              140.0);
    	streetTypes.put("motorway_link",          50.0);
        streetTypes.put("path",                    4.5);
        streetTypes.put("pedestrian",              4.5);
        streetTypes.put("primary",               100.0);
        streetTypes.put("primary_link",           50.0);
        streetTypes.put("residential",            30.0);
        streetTypes.put("road",                   70.0);
        streetTypes.put("secondary",             100.0);
        streetTypes.put("secondary_link",         50.0);
        streetTypes.put("service",                50.0);
        streetTypes.put("steps",                   3.0);
        streetTypes.put("tertiary",               70.0);
        streetTypes.put("tertiary_link",          50.0);
        streetTypes.put("track",                  50.0);
        streetTypes.put("track:grade1",           50.0);
        streetTypes.put("track:grade2",           30.0);
        streetTypes.put("track:grade3",           20.0);
        streetTypes.put("trunk",                 120.0);
        streetTypes.put("trunk_link",             50.0);
        streetTypes.put("unclassified",           50.0);
    }

    private Util() {}

    /**
     * Writes an array of strings to a file. Every entry of the array starts in a new line.
     *
     * @param s             Array of strings
     * @param outputFile    The file to write in
     */
    public static void saveStrings(String[] s, File outputFile) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            for (String t : s) {
                writer.write(t + "\n");
            }
            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println("Error while writing");
            e.printStackTrace();
        }
    }

    /**
     * Writes a graph starting with its number of nodes and edges, followed by a list of all nodes and finally a list
     * of all edges.
     *
     * @param nodes             Array of nodes to print
     * @param edges             List of edges to print
     * @param graphFile         The file to write in
     * @param coordinateFile    The file to write in
     */
    public static void saveGraph(Node[] nodes, List<Edge> edges, File graphFile, File coordinateFile) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(graphFile));
            writer.write("p sp " + nodes.length + " " + edges.size());
            writer.newLine();
            int i = 0;
            for (Edge edge : edges) {
                writer.write(edge.toString());
                writer.newLine();
                if (i % 1000 == 0) {
                   writer.flush();
                }
                i++;
            }
            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println("Error while writing");
            e.printStackTrace();
        }
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(coordinateFile));
            writer.write("p aux sp co " + nodes.length);
            writer.newLine();
            for (Node node : nodes) {
                writer.write(node.toString());
                writer.newLine();
                if (node.getGraphId() % 1000 == 0) {
                   writer.flush();
                }
            }
            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println("Error while writing");
            e.printStackTrace();
        }
    }

    /**
     * Writes a list of nodes.
     *
     * @param nodes       List of nodes to print
     * @param nodeFile    The file to write in
     */
    public static void saveNodes(List<Node> nodes, File nodeFile) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(nodeFile));
            for (Node node : nodes) {
                writer.write(node.toString());
                writer.newLine();
            }
            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println("Error while writing");
            e.printStackTrace();
        }
    }

    /**
     * Opens a file dialog, so that the user can specify a file to read or write in.
     *
     * @param fileSufix     The file-type the user should select
     * @param mode          Specifies if the file should be opened or saved
     * @return              The file specified by the user
     */
    public static File askFile(final String fileSufix, int mode) {
        System.out.println("\tasking for file");
        JFileChooser fc = new JFileChooser();
        fc.setFileFilter(new FileFilter() {
            @Override
            public String getDescription() {
                return fileSufix;
            }

            @Override
            public boolean accept(File f) {
                return f.isDirectory() || f.getName().toLowerCase().matches("^.*\\" + fileSufix + "$");
            }
        });
        fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        if (mode == OPEN) {
            fc.showOpenDialog(null);
        } else if (mode == SAVE) {
            fc.showSaveDialog(null);
        }
        File result = fc.getSelectedFile();
        System.out.println("\trecived file");
        return result;
    }
}
