/* $RCSfile$
 * $Author: hansonr $
 * $Date: 2006-08-02 11:48:43 -0500 (Wed, 02 Aug 2006) $
 * $Revision: 5364 $
 *
 * Copyright (C) 2003-2005  Miguel, Jmol Development, www.jmol.org
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package com.actelion.research.chem.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.Hashtable;
import java.util.Map;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXParseException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Adapted from org.jmol.adapters.readers.xml.XmlReader.
 * 
 * A generic XML reader template -- by itself, does nothing.
 * 
 * XmlReader takes all XML streams, whether from a file reader or from DOM.
 * 
 * This class functions as a resolver, since it: (1) identifying the specific
 * strain of XML to be handled, and (2) passing the responsibility on to the
 * correct format-specific XML readers. There are parallel entry points and
 * handler methods for reader and DOM. Each format-specific XML reader then
 * assigns its own handler to manage the parsing of elements.
 * 
 * In addition, this class handles generic XML tag parsing.
 * 
 * XmlHandler extends DefaultHandler is the generic interface to both reader and
 * DOM element parsing.
 * 
 * Note that the tag processing routines are shared between SAX and DOM
 * processors. This means that attributes must be transformed from either
 * Attributes (SAX) or JSObjects (DOM) to Hashtable name:value pairs. This is
 * taken care of in XmlHandler for all readers.
 * 
 * Test files:
 * 
 * molpro: vib.xml odyssey: water.xodydata cml: a wide variety of files in
 * data-files.
 * 
 * -Bob Hanson 2006.03.27
 * 
 */

abstract public class XmlReader {

  public Map<String, String> atts = new Hashtable<String, String>();

  protected String err;
  	
  protected BufferedReader reader;

  protected String parseXML() {
    org.xml.sax.XMLReader saxReader = null;
    /**
     * @j2sNative
     * 
     * 
     * 
     */
    {
      try {
        javax.xml.parsers.SAXParserFactory spf = javax.xml.parsers.SAXParserFactory
            .newInstance();
        spf.setNamespaceAware(true);
        javax.xml.parsers.SAXParser saxParser = spf.newSAXParser();
        saxReader = saxParser.getXMLReader();

        // otherwise, DTD UTI is treated as a URL, retrieved, and scanned.
        // see https://stackoverflow.com/questions/10257576/how-to-ignore-inline-dtd-when-parsing-xml-file-in-java
        saxReader.setFeature("http://xml.org/sax/features/validation", false);
        saxReader.setFeature("http://apache.org/xml/features/nonvalidating/load-dtd-grammar", false);
        saxReader.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
        
        System.out.println("Using JAXP/SAX XML parser.");
      } catch (Exception e) {
          System.out.println("Could not instantiate JAXP/SAX XML reader: " + e.getMessage());
      }
      if (saxReader == null)
        return "No XML reader found";
    }
	try {
		processXml(saxReader);
	    return null;
	} catch (Exception e) {
		e.printStackTrace();
		return "Error reading XML: " + (e.getMessage());

	}
  }

  /**
   * @param saxReader
   * @throws Exception
   */
  protected void processXml(Object saxReader)
      throws Exception {
    BufferedReader rdr = previewXML(reader);
    if (saxReader == null) {
    	parseJSDOM(rdr);
    } else {
    	new XmlHandler().parseXML(this, saxReader, rdr);
    }
  }

  private void parseJSDOM(BufferedReader rdr) throws UnsupportedEncodingException {
      attribs = new Object[1];
      domObj = new Object[1];
      Object o = "";
      byte[] data = null;
      /**
       * 
       * @j2sNative
       * 
       *            o = rdr.lock.lock; 
       *            if (o && o.$in) {
       *              data = o.$in.buf;
       *            } else if ((o=rdr.$in.str) != null) {
       *            } else if (rdr.$in.$in.$in.fd) {
       *              // may need to adjust this;
       *              o = rdr.$in.$in;
       *              data = o.$in.fd._file.秘bytes; 
       *            } else {
       *              data = (o=rdr.$in.$in).$in.buf;
       *            }
       */
      {
      }
      if (o instanceof BufferedInputStream)
        o = new String(data, "utf-8");
      boolean isjs = false;
      /**
       * 
       * @j2sNative
       * 
       *            isjs = true;
       * 
       */
      {
        walkDOMTree();
      }
      if (isjs) {
        this.domObj[0] = this.createDomNodeJS("xmlReader", o);
        this.walkDOMTree();
        this.createDomNodeJS("xmlReader", null);
      }
  }

/**
   * An opportunity to fix XML (nmrML issues)
   * @param reader
   * @return reader or a repaired reader
   * @throws IOException 
   */
  protected BufferedReader previewXML(BufferedReader reader) throws IOException {
    return reader;
  }

  /**
   * @param id  
   * @param data 
   * @return dom object 
   */
  Object createDomNodeJS(String id, Object data) {
    // no doubt there is a more efficient way to do this.
    // Firefox, at least, does not recognize "/>" in HTML blocks
    // that are added this way.
    
    Object d = null;
    /**
     * note that there is no need to actually load it into the document
     * 
     * @j2sNative
     * 
      if (!data)
        return null;
      if (data.indexOf("<?") == 0)
        data = data.substring(data.indexOf("<", 1));
      if (data.indexOf("/>") >= 0) {
        var D = data.split("/>");
        for (var i = D.length - 1; --i >= 0;) {
          var s = D[i];
          var pt = s.lastIndexOf("<") + 1;
          var pt2 = pt;
          var len = s.length;
          var name = "";
          while (++pt2 < len) {
            if (" \t\n\r".indexOf(s.charAt(pt2))>= 0) {
              var name = s.substring(pt, pt2);
              D[i] = s + "></"+name+">";
              break;
            }     
          }
        }
        data = D.join('');
      }
      d = document.createElement("_xml");
      d.innerHTML = data;
     * 
     */
    {
      // only called by j2s
    }
    return d;
  }
  
  ////////////////////////////////////////////////////////////////

  /**
   * 
   * @param localName
   * @param nodeName
   */
  abstract protected void processStartElement(String localName, String nodeName);
  /*
   *  keepChars is used to signal 
   *  that characters between end tags should be kept
   *  
   */

  protected boolean keepChars;
  protected StringBuffer chars = new StringBuffer();
  
  protected void setKeepChars(boolean TF) {
    keepChars = TF;
    chars.setLength(0);
  }

  /**
   * 
   * @param localName
   */
  abstract void processEndElement(String localName);

  //////////////////// DOM or JavaScript parsing /////////////////

  // walk DOM tree given by JSObject. For every element, call
  // startElement with the appropriate strings etc., and then
  // endElement when the element is closed.

  private Object[] domObj = new Object[1];
  private Object[] attribs;
//  private Object[] attArgs;
//  private Object[] nullObj = new Object[0];

  private void walkDOMTree() {
    String localName;
    /**
     * @j2sNative
     * 
     * localName = "nodeName";
     * 
     */
    {
      localName = "localName";
    }
    String nodeName = ((String) jsObjectGetMember(domObj, localName));
    localName = fixLocal(nodeName);
    if (localName == null)
      return;
    if (localName.equals("#text")) {
      if (keepChars)
        chars.append((String) jsObjectGetMember(domObj, "data"));
      return;
    }
    localName = localName.toLowerCase();
    nodeName = nodeName.toLowerCase();
    attribs[0] = jsObjectGetMember(domObj, "attributes");
    getDOMAttributesA(attribs);
    processStartElement(localName, nodeName);
    boolean haveChildren = false;
    /**
     * @j2sNative
     * 
     *            haveChildren = this.domObj[0].hasChildNodes;
     * 
     */
    {
//      haveChildren = ((Boolean) jsObjectCall(domObj, "hasChildNodes", null))
//          .booleanValue();
    }
    if (haveChildren) {
      Object nextNode = jsObjectGetMember(domObj, "firstChild");
      while (nextNode != null) {
        domObj[0] = nextNode;
        walkDOMTree();
        domObj[0] = nextNode;
        nextNode = jsObjectGetMember(domObj, "nextSibling");
      }
    }
    processEndElement(localName);
  }

  private String fixLocal(String name) {
    /**
     * @j2sNative
     * 
     *            var pt = (name== null ? -1 : name.indexOf(":")); return (pt >=
     *            0 ? name.substring(pt+1) : name);
     */
    {
      return name;
    }
  }

  private class NVPair {
    String name;
    String value;
  }
  
  @SuppressWarnings("unused")
  private void getDOMAttributesA(Object[] attributes) {

    atts.clear();
    if (attributes == null)
      return;

    NVPair[] nodes = null;

    /**
     * @j2sNative
     * 
     * 
     *            nodes = attributes[0];
     * 
     */
    {
      if (true)
        return;
    }
    // JavaScript only
    for (int i = nodes.length; --i >= 0;)
      atts.put(fixLocal(nodes[i].name).toLowerCase(), nodes[i].value);

  }
  
  

  /**
   * @param jsObject  
   * @param name 
   * @return an object
   */
  private Object jsObjectGetMember(Object[] jsObject, String name) {
    /**
     * @j2sNative
     * 
     * return jsObject[0][name]; 
     * 
     */
    {
    return null;
    }
  }

  public void endDocument() {
	  // can be implmented
  }

  public static class XmlHandler extends DefaultHandler {

	  private XmlReader xmlReader;

	  public XmlHandler() {
	    // for reflection
	  }
	  
	  void parseXML(XmlReader xmlReader, Object saxReaderObj, BufferedReader reader) throws Exception {
	    this.xmlReader = xmlReader;
	    XMLReader saxReader = (XMLReader) saxReaderObj;
	    saxReader.setFeature("http://xml.org/sax/features/validation", false);
	    saxReader.setFeature("http://xml.org/sax/features/namespaces", true);
	    saxReader.setEntityResolver(this);
	    saxReader.setContentHandler(this);
	    saxReader.setErrorHandler(this);
	    InputSource is = new InputSource(reader);
	    is.setSystemId("foo");
	    saxReader.parse(is);
	  }

	  @Override
	  public void startDocument() {
	  }

	  @Override
	  public void endDocument() {
	    xmlReader.endDocument();
	  }

	//  /*
	//   * see org.xml.sax.ContentHandler#startElement(java.lang.String, java.lang.String, java.lang.String, org.xml.sax.Attributes)
	//   * startElement and endElement should be extended in each reader
	//   */
	//
	//  private String debugContext = "";

	  @Override
	  public void startElement(String namespaceURI, String localName, String nodeName,
	                           Attributes attributes) {
	    // Note added 5/2015 BH 
	    //
	    // There seems to be an issue with what is the localName and what is the qName.
	    // For example, for
	    //
	    //    <cml:molecule xmlns:cml="http://www.xml-cml.org/schema/cml2/core">
	    //
	    // Both HTML5 xmlDocument and JAXPSAXParser report:
	    //
	    //  namespaceURI = "http://www.xml-cml.org/schema/cml2/core"
	    //  localName = "molecule"
	    //  nodeName = "cml:molecule
	    //
	    //  But com.sun.apache.xerces.internal.jaxp.SAXParserImpl$JAXPSAXParser reports:
	    //
	    //  namespaceURI = ""
	    //  localName = ""
	    //  nodeName = "cml:molecule"
	    //
	    // You would think that is a bug...
	    //
	    // I only realized this recently when I wrote the JavaScript swingjs.JSSAXParser
	    // code and tested it against javax.xml.parsers.SAXParserFactory.newInstance().newSAXParser().
	    //
	    xmlReader.atts.clear();
	    for (int i = attributes.getLength(); --i >= 0;)
	      xmlReader.atts.put(attributes.getLocalName(i).toLowerCase(), attributes.getValue(i));
	    xmlReader.processStartElement(localName.toLowerCase(), nodeName.toLowerCase());
	  }

	  @Override
	  public void endElement(String uri, String localName, String qName) {
	    xmlReader.processEndElement(localName.toLowerCase());
	  }

	  @Override
	  public void characters(char[] ch, int start, int length) {
	    if (xmlReader.keepChars)
	      xmlReader.chars.append(ch, start, length);
	  }

	  // Methods for entity resolving, e.g. getting a DTD resolved

	  public InputSource resolveEntity(String name, String publicId,
	                                   String baseURI, String systemId) {
	    return null;
	  }

	  @Override
	  public void error(SAXParseException exception) {
	    System.err.println("SAX ERROR:" + exception.getMessage());
	  }
	  @Override
	  public void fatalError(SAXParseException exception) {
	    System.err.println("SAX FATAL:" + exception.getMessage());
	  }

	  @Override
	  public void warning(SAXParseException exception) {
		  System.err.println("SAX WARNING:" + exception.getMessage());
	  }

	}
  
}
