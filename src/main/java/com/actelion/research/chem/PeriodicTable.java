/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.chem;

import java.util.Hashtable;

public class PeriodicTable {

	public static final int ConnectionPoint = 0;

	public static final int Hydrogen = 1;

	public static final int Carbon = 6;
	
	public static final int Nitrogen = 7;
	
	public static final int Oxygen = 8;
	
	public static final int Fluorine = 9;
	
	public static final int Silicon = 14;
	
	public static final int Phosphorus = 15;
	
	public static final int Sulfur = 16;
	
	public static final int Chlorine = 17;
	
	public static final int Bromine = 35;
	
	public static final int Iodine = 53;
	
	
	
	private static Element[] arrData;

	private Hashtable<Integer, Element> htblDataAtNo;

	private Hashtable<String, Element> htblDataName;

	private Hashtable<String, Element> htblDataSymbol;

	private static PeriodicTable tbl;

	private PeriodicTable() {

		htblDataAtNo = new Hashtable<Integer, Element>();

		htblDataName = new Hashtable<String, Element>();

		htblDataSymbol = new Hashtable<String, Element>();

// Order Number, Element Symbol, Atomic Weight, Covalent Radius, VdW radius
    arrData = new Element [] {
        new Element(1, "Hydrogen"      , "H"  ,1.008	,0.435	,1.185	),
        new Element(2, "Helium"        , "He" ,4.003	,0.93	,1.785	),
        new Element(3, "Lithium"       , "Li" ,6.94	    ,1.5199	,0	),
        new Element(4, "Beryllium"     , "Be" ,9.013	,1.143	,0	),
        new Element(5, "Boron"         , "B"  ,10.82	,0.975	,1.7	),
        new Element(6, "Carbon"        , "C"  ,12.011	,0.655	,1.75	),
        new Element(7, "Nitrogen"      , "N"  ,14.008	,0.75	,1.525	),
        new Element(8, "Oxygen"        , "O"  ,16	    ,0.73	,1.4	),
        new Element(9, "Fluorine"      , "F"  ,19	    ,0.72	,1.35	),
        new Element(10, "Neon"         , "Ne" ,20.183	,0.71	,1.6	),
        new Element(11, "Sodium"       , "Na" ,22.991	,1.8579	,0	),
        new Element(12, "Magnesium"    , "Mg" ,24.32	,1.6047	,0	),
        new Element(13, "Aluminum"     , "Al" ,26.98	,1.4318	,0	),
        new Element(14, "Silicon"      , "Si" ,28.09	,1.1758	,2	),
        new Element(15, "Phosphorus"   , "P"  ,30.975	,1.06	,1.9	),
        new Element(16, "Sulfur"       , "S"  ,32.066	,1.02	,1.85	),
        new Element(17, "Chlorine"     , "Cl" ,35.457	,0.99	,1.78	),
        new Element(18, "Argon"        , "Ar" ,39.944	,0.98	,1.9	),
        new Element(19, "Potassium"    , "K"  ,39.1	    ,2.262	,0	),
        new Element(20, "Calcium"      , "Ca" ,40.08	,1.9758	,0	),
        new Element(21, "Scandium"     , "Sc" ,44.96	,1.6545	,0	),
        new Element(22, "Titanium"     , "Ti" ,47.9	    ,1.4755	,0	),
        new Element(23, "Vanadium"     , "V"  ,50.95	,1.309	,0	),
        new Element(24, "Chromium"     , "Cr" ,52.01	,1.249	,0	),
        new Element(25, "Manganese"    , "Mn" ,54.94	,1.35	,0	),
        new Element(26, "Iron"         , "Fe" ,55.85	,1.2411	,0	),
        new Element(27, "Cobalt"       , "Co" ,58.94	,1.2535	,0	),
        new Element(28, "Nickel"       , "Ni" ,58.69	,1.246	,0	),
        new Element(29, "Copper"       , "Cu" ,63.54	,1.278	,0	),
        new Element(30, "Zinc"         , "Zn" ,65.38	,1.3325	,0	),
        new Element(31, "Gallium"      , "Ga" ,69.72	,1.3501	,0	),
        new Element(32, "Germanium"    , "Ge" ,72.6	    ,1.2248	,0	),
        new Element(33, "Arsenic"      , "As" ,74.91	,1.2	,2.1	),
        new Element(34, "Selenium"     , "Se" ,78.96	,1.16	,2	),
        new Element(35, "Bromine"      , "Br" ,79.916	,1.14	,1.97	),
        new Element(36, "Krypton"      , "Kr" ,83.8	    ,1.12	,2	),
        new Element(37, "Rubidium"     , "Rb" ,85.48	,2.47	,0	),
        new Element(38, "Strontium"    , "Sr" ,87.63	,2.1513	,0	),
        new Element(39, "Yttrium"      , "Y"  ,88.92	,1.8237	,0	),
        new Element(40, "Zirconium"    , "Zr" ,91.22	,1.6156	,0	),
        new Element(41, "Niobium"      , "Nb" ,92.91	,1.4318	,0	),
        new Element(42, "Molybdenum"   , "Mo" ,95.95	,1.3626	,0	),
        new Element(43, "Technetium"   , "Tc" ,99	    ,1.3675	,0	),
        new Element(44, "Ruthenium"    , "Ru" ,101.1	,1.3529	,0	),
        new Element(45, "Rhodium"      , "Rh" ,102.91	,1.345	,0	),
        new Element(46, "Palladium"    , "Pd" ,106.7	,1.3755	,0	),
        new Element(47, "Silver"       , "Ag" ,107.88	,1.4447	,0	),
        new Element(48, "Cadmium"      , "Cd" ,112.41	,1.4894	,0	),
        new Element(49, "Indium"       , "In" ,114.76	,1.6662	,0	),
        new Element(50, "Tin"          , "Sn" ,118.7	,1.5375	,0	),
        new Element(51, "Antimony"     , "Sb" ,121.76	,1.4	,2.2	),
        new Element(52, "Tellurium"    , "Te" ,127.61	,1.36	,2.2	),
        new Element(53, "Iodine"       , "I"  ,126.91	,1.33	,2.2	),
        new Element(54, "Xenon"        , "Xe" ,131.3	,1.31	,2.2	),
        new Element(55, "Cesium"       , "Cs" ,132.91	,2.6325	,0	),
        new Element(56, "Barium"       , "Ba" ,137.36	,2.1705	,0	),
        new Element(57, "Lanthanum"    , "La" ,138.92	,1.8725	,0	),
        new Element(58, "Cerium"       , "Ce" ,140.13	,1.8243	,0	),
        new Element(59, "Praesodymium" , "Pr" ,140.92	,1.8362	,0	),
        new Element(60, "Neodymium"    , "Nd" ,144.27	,1.8295	,0	),
        new Element(61, "Promethium"   , "Pm" ,145	    ,1.809	,0	),
        new Element(62, "Samarium"     , "Sm" ,150.43	,1.804	,0	),
        new Element(63, "Europium"     , "Eu" ,152	    ,1.984	,0	),
        new Element(64, "Gadolinium"   , "Gd" ,156.9	,1.818	,0	),
        new Element(65, "Terbium"      , "Tb" ,158.93	,1.8005	,0	),
        new Element(66, "Dyprosium"    , "Dy" ,162.46	,1.7951	,0	),
        new Element(67, "Holmium"      , "Ho" ,164.94	,1.7886	,0	),
        new Element(68, "Erbium"       , "Er" ,167.2	,1.7794	,0	),
        new Element(69, "Thulium"      , "Tm" ,168.94	,1.7687	,0	),
        new Element(70, "Ytterbium"    , "Yb" ,173.04	,1.9396	,0	),
        new Element(71, "Lutetium"     , "Lu" ,174.99	,1.7515	,0	),
        new Element(72, "Hafnium"      , "Hf" ,178.6	,1.5973	,0	),
        new Element(73, "Tantalium"    , "Ta" ,180.95	,1.428	,0	),
        new Element(74, "Wolfram"      , "W"  ,183.92	,1.3705	,0	),
        new Element(75, "Rhenium"      , "Re" ,186.31	,1.38	,0	),
        new Element(76, "Osmium"       , "Os" ,190.2	,1.3676	,0	),
        new Element(77, "Iridium"      , "Ir" ,192.2	,1.3573	,0	),
        new Element(78, "Platinum"     , "Pt" ,195.23	,1.3873	,0	),
        new Element(79, "Gold"         , "Au" ,197	    ,1.4419	,0	),
        new Element(80, "Mercury"      , "Hg" ,200.61	,1.5025	,0	),
        new Element(81, "Thallium"     , "Tl" ,204.39	,1.7283	,0	),
        new Element(82, "Lead"         , "Pb" ,207.21	,1.7501	,0	),
        new Element(83, "Bismuth"      , "Bi" ,209	    ,1.46	,0	),
        new Element(84, "Polonium"     , "Po" ,210	    ,1.46	,0	),
        new Element(85, "Astatine"     , "At" ,210	    ,1.45	,0	),
        new Element(86, "Radon"        , "Rn" ,222	    ,1.43	,0	),
        new Element(87, "Francium"     , "Fr" ,223	    ,2.5	,0	),
        new Element(88, "Radium"       , "Ra" ,226.05	,2.14	,0	),
        new Element(89, "Actinium"     , "Ac" ,227	    ,1.8775	,0	),
        new Element(90, "Thorium"      , "Th" ,232.05	,1.7975	,0	),
        new Element(91, "Protactinium" , "Pa" ,231	    ,1.6086	,0	),
        new Element(92, "Uranium"      , "U"  ,238.07	,1.5683	,0	),
        new Element(93, "Neptunium"    , "Np" ,237	    ,0	    ,0	),
        new Element(94, "Plutonium"    , "Pu" ,242	,0	,0	),
        new Element(95, "Americium"    , "Am" ,243	,0	,0	),
        new Element(96, "Curium"       , "Cm" ,243	,0	,0	),
        new Element(97, "Berkelium"    , "Bk" ,245	,0	,0	),
        new Element(98, "Californium"  , "Cf" ,246	,0	,0	),
        new Element(99, "Einsteinium"  , "E"  ,-999	,0	,0	),
        new Element(100, "Fermium"     , "Fm" ,-999	,0	,0	),
        new Element(101, "Mendelevium" , "Mv" ,-999	,0	,0	),
        new Element(120, "Dummy"       , "*"  ,0.0	,0  ,0      ),
        new Element(121, "Any"         , "A"  ,0.0	,0	,0      ),
        new Element(122, "List"        , "L"  ,0.0	,0	,0      )};

    for (int ii = 0; ii < arrData.length; ii++) {
			htblDataAtNo
					.put(new Integer(arrData[ii].getOrderNumber()), arrData[ii]);
			htblDataName.put(arrData[ii].getName(), arrData[ii]);
			htblDataSymbol.put(arrData[ii].getSymbol(), arrData[ii]);
		}
	}

	private static PeriodicTable getInstance() {
		if (tbl == null) {
			tbl = new PeriodicTable();
		}
		return tbl;
	}

	/**
	 * 
	 * @param sNameOrSymbol
	 * @return atomic order number or 0 if the input can not be resolved to an element.
	 */
	public static int number(String sNameOrSymbol) {

		int iOrderNumber = 0;

		Element el = getInstance().htblDataName.get(sNameOrSymbol);
		if (el == null)
			el = getInstance().htblDataSymbol.get(sNameOrSymbol);

		if (el != null)
			iOrderNumber = el.getOrderNumber();

		return iOrderNumber;
	}

	public static String symbol(int iOrderNumber) {

		String sSymbol = "";

		Element el = getInstance().htblDataAtNo.get(iOrderNumber);

		if (el != null)
			sSymbol = el.getSymbol();

		return sSymbol;
	}

	public static int size(){
		return getInstance().htblDataAtNo.size();
	}
	
	/**
	 * 
	 * @param iOrderNumber
	 * @return element weight
	 */
	public static double getWeight(int iOrderNumber) {

		Element el = getInstance().htblDataAtNo.get(iOrderNumber);

		return el.getWeight();
	}

	public static String name(int iOrderNumber) {

		String sName = "";

		Element el = getInstance().htblDataAtNo.get(iOrderNumber);

		if (el != null)
			sName = el.getName();

		return sName;
	}

	public static Element getElement(int orderNumber) {
		return getInstance().htblDataAtNo.get(orderNumber);
	}
	
	public String toString() {
		String sPeriodicTable = "";
		return sPeriodicTable;
	}
}