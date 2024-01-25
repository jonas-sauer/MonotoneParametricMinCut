package OSM;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import static javax.xml.stream.XMLStreamReader.*;

public class XMLParser {
    
    public static XMLInputFactory factory = XMLInputFactory.newInstance();
    
    private Map<String, XMLNodeListener> listeners = new HashMap<String, XMLParser.XMLNodeListener>();
    
    public void addNodeListener(XMLNodeListener l) {
        this.listeners.put(l.getNodeName(), l);
    }
    
    public void removeNodeListener(XMLNodeListener l ) {
        this.listeners.remove(l.getNodeName());
    }
    
    public void process(File xml) {
        try {
            FileReader fr = new FileReader(xml);
            XMLStreamReader xmlReader = factory.createXMLStreamReader(fr);
            processTag(xmlReader, this.listeners);
            xmlReader.close();
            fr.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (XMLStreamException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private void processTag(XMLStreamReader xmlReader, Map<String, XMLNodeListener> listeners) throws XMLStreamException {
        while (xmlReader.hasNext()) {
            int type = xmlReader.next();
            if (type == START_ELEMENT) {
                String tagName = xmlReader.getName().toString();
                if (listeners.containsKey(tagName)) {
                    XMLNodeListener l = listeners.get(tagName);
                    l.nodeBegin(new Element(xmlReader));
                    processTag(xmlReader, l.getSubNodeListeners());
                    l.nodeEnd();
                } else {
                    processTag(xmlReader, listeners);
                }
            } else if (type == END_ELEMENT) {
                return;           
            }
        }
    }

    public static interface XMLNodeListener {
        public String getNodeName();
        public Map<String, XMLNodeListener> getSubNodeListeners();
        public void nodeBegin(Element e);
        public void nodeEnd();
    }
    
    public static abstract class InnerNodeListener implements XMLNodeListener {
        private Map<String, XMLNodeListener> subNodeListeners = new HashMap<String, XMLNodeListener>();
        private String tag;
        public InnerNodeListener(String tag) {
            this.tag = tag;
        }
        public Map<String, XMLNodeListener> getSubNodeListeners() {
            return subNodeListeners;
        }
        public String getNodeName() {
            return tag;
        }
        public void nodeEnd() {}
    }
    
    public static final class Element {
        private XMLStreamReader reader;
        public Element(XMLStreamReader reader) {
            this.reader = reader;
        }   
        public boolean hasAttribute(String name, String value) {
            return value.equals(getAttribute(name, null));
        }
        public String getAttribute(String key) {
            return getAttribute(key, "");
        }        
        public String getAttribute(String key, String defaultValue) {
            for (int i = 0; i < reader.getAttributeCount(); i++) {
                if (reader.getAttributeName(i).toString().equals(key)) {
                    defaultValue = reader.getAttributeValue(i);
                }
            }
            return defaultValue;
        }        
    }
    
}
