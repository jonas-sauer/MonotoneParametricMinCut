package OSM;

import java.io.File;
import java.io.IOException;

import static OSM.Parser.*;
import static OSM.Util.*;

public class Main {

    public static void main(String[] args) {
    	CommandLineParser clp = new CommandLineParser(args);
        if (clp.getBoolean("h")) help();

        System.out.println("Start");
        File input = new File(clp.getString("i", ""));
        if (!input.exists() || !input.canRead()) {
            input = askFile(".osm", OPEN);
        }
        if (input != null) {
            System.out.println("\tinput file is: " + input.getAbsolutePath());
            Graph graph = parseOSM(input);
            File outputGr = new File(clp.getString("gr", ""));
            File outputCo = new File(clp.getString("co", ""));
            if (!outputGr.exists()) {
                try {
                	outputGr.createNewFile();
                } catch (IOException e) {
                    System.out.println("\tERROR! Could not create graph file: " + outputGr.getAbsolutePath());
                }
            }
            if (!outputCo.exists()) {
                try {
                	outputCo.createNewFile();
                } catch (IOException e) {
                    System.out.println("\tERROR! Could not create coordinates file: " + outputCo.getAbsolutePath());
                }
            }
            if (!outputGr.exists() || !outputGr.canWrite()) {
            	outputGr = askFile(".gr", SAVE);
            }
            if (!outputCo.exists() || !outputCo.canWrite()) {
            	outputCo = askFile(".co", SAVE);
            }
            if (outputGr != null && outputCo != null) {
                System.out.println("\toutput file is: " + outputGr.getAbsolutePath() + ", " + outputCo.getAbsolutePath());
                saveGraph(graph.nodes, graph.edges, outputGr, outputCo);
                System.out.println("Done.");
            } else {
                System.out.println("No output file");
            }
        } else {
            System.out.println("No input file");
        }
    }

    public static void help() {
        System.out.println("Usage: java " + Main.class.getName() + "<options>");
        System.out.println("-h\t\t -- show help");
        System.out.println("-i <filename>\t -- input file (.osm file)");
        System.out.println("-gr <filename>\t -- output file (.gr file)");
        System.out.println("-co <filename>\t -- output file (.co file)");
        System.exit(0);
    }

}
