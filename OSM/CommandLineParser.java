package OSM;

public class CommandLineParser {
    
    private String[] args;

    public CommandLineParser(String[] args) {
        this.args = args;
    }
    
    public String getString(String key) {
        return getString(key, null);
    }

    public String getString(String key, String defaultValue) {
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].equals("-" + key)) {
                return args[i + 1];
            }
        }
        return defaultValue;
    }
    
    public boolean getBoolean(String key) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-" + key)) {
                return true;
            }
        }
        return false;
    }
    
    public int getInt(String key) {
        return Integer.parseInt(getString(key, "0"));
    }
    
    public int getInt(String key, int defaultValue) {
        return Integer.parseInt(getString(key, "" + defaultValue));
    }
    
    public double getDouble(String key) {
        return Double.parseDouble(getString(key, "0"));
    }
    
    public double getDouble(String key, double defaultValue) {
        return Double.parseDouble(getString(key, "" + defaultValue));
    }
        
}
