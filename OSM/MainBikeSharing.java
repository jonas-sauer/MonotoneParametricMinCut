package OSM;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static OSM.Parser.*;
import static OSM.Util.*;

public class MainBikeSharing {

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
            List<Node> stations = parseBickeSharing(input);
            File outputTxt = new File(clp.getString("o", ""));
            if (!outputTxt.exists()) {
                try {
                    outputTxt.createNewFile();
                } catch (IOException e) {
                    System.out.println("\tERROR! Could not create graph file: " + outputTxt.getAbsolutePath());
                }
            }
            if (!outputTxt.exists() || !outputTxt.canWrite()) {
                outputTxt = askFile(".txt", SAVE);
            }
            if (outputTxt != null) {
                System.out.println("\toutput file is: " + outputTxt.getAbsolutePath());
                saveNodes(stations, outputTxt);
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
        System.out.println("-o <filename>\t -- output file (.txt file)");
        System.exit(0);
    }

}
