package pepXML2OUT_jmb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.*;                   
import java.util.Scanner;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Vector;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.XMLReader;

import pepXML2OUT_jmb.DocHandlers;

public class pepXML2OUT_jmb {
    private NodeList nodelist;
    private String basef="", fname="";
    String XmlFile;
    
    public static void main(String[] args) throws IOException, ParserConfigurationException, SAXException{ // added IOException
        if(args.length!=1){
            System.out.println("Error: Invalid number of parameter(s)");
            return;
        } 
        System.out.println(args[0]);
        
        pepXML2OUT_jmb p  = new pepXML2OUT_jmb();
        p.setFileNames(args[0]);
        
        DocHandlers handler = new DocHandlers();
        handler.setBaseFile(p.getBaseF());
        handler.setFileName(p.getFName());
        
        try {
            SAXParserFactory spf = SAXParserFactory.newInstance();
            spf.setNamespaceAware(true);
            SAXParser saxParser = spf.newSAXParser();
            saxParser.parse( convertToFileURL(args[0]) ,handler);
        } catch (ParserConfigurationException e) {
            System.out.println("ParserConfig error");
        } catch (SAXException e) {
            System.out.println("SAXException : xml not well formed");
        }
        
        //XMLReader xmlReader = saxParser.getXMLReader();
        //xmlReader.setContentHandler(new DocHandlers());
                
        //xmlReader.parse(convertToFileURL(args[0]));
        /* Change the first line of input file from "UTF-8" to "ISO8859-1" before p.loadXML() */
        
    }       
        
private static String convertToFileURL(String filename) {
        String path = new File(filename).getAbsolutePath();
        if (File.separatorChar != '/') {
            path = path.replace(File.separatorChar, '/');
        }

        if (!path.startsWith("/")) {
            path = "/" + path;
        }
        return "file:" + path;
    }   
void setFileNames(String file){
    basef=file.substring(0,file.lastIndexOf("/")>0?file.lastIndexOf("/")+1:file.lastIndexOf("\\")+1);
    fname=file.substring(file.lastIndexOf("/")>0?file.lastIndexOf("/")+1:file.lastIndexOf("\\")+1,file.lastIndexOf("."));
    }
String getFName(){
        return fname;
    }
String getBaseF(){
        return basef;
    }
}        
        