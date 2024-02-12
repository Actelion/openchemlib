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

package com.actelion.research.chem.forcefield.mmff;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;

/**
 * Basic CSV parser. The Csv class provides very basic CSV file parsing
 * and assumes that CSV files are in a strict format. CSV files are
 * expected to conform to the following syntax:
 *
 *  - First line contains only the number of rows to be read as an
 *      integer.
 *
 *  - Second line is a strictly comma separated row of type strings, where
 *      each type string has no preceeding or following whitespace. An
 *      example of a type string would be:
 *
 *          int,int,float,char
 *
 *      This type string informs the parser that the rows of the table
 *      have four cells of which the first two are integers and the last
 *      two are floating point numbers. Only three type strings are
 *      supported: 'int', 'float' and 'char'.
 *
 *  - Third line and onwards. These are the table rows. They are expected
 *      to conform to the type format specified in the type line (line 2).
 *      As well as this all rows are expected to have length equal to the
 *      type row. Cells in each row are strictly comma separated and
 *      contain no whitespace. Each cell is expected to be in a format
 *      that Java can parse using parseInt() for integer types,
 *      parseDouble() for floating point types, and char data types are
 *      expected to be a single character optionally enclosed in single
 *      or double quotation marks. An example of several CSV rows which
 *      adhere to the previously given type row:
 *
 *          1,16,3.900,'A'
 *          1,17,2.700,'-'
 *          1,18,2.100,'B'
 *          2,19,4.500,'C'
 *          2,20,1.050,'-'
 *          2,21,0.150,'D'
 */
public final class Csv {
    /**
     * Read a file at location 'path' and parse it as a CSV file.
     *  @param path The string path to the CSV file.
     *  @return A two dimensional object array where each row is one row
     *      from the CSV file and each sub array element is one cell from
     *      one row in the CSV file. The underlying objects are of the
     *      correct boxed type as specified in the type line of the CSV
     *      file and can be cast back to their type.
     */
    public static Object[][] readFile(String path) {
        BufferedReader br = null;

        try {
            String line;
            br = new BufferedReader(new InputStreamReader(Csv.class.getResourceAsStream(path), StandardCharsets.UTF_8));

            int size = Integer.parseInt(br.readLine().trim());
            String[] format = br.readLine().trim().split(",");
            Object[][] table = new Object[size][format.length];

            for (int ln=0; (line = br.readLine()) != null && ln<size; ln++) {
                String[] row = line.trim().split(",");

                for (int c=0; c<format.length; c++) {
                    switch (format[c].charAt(0)) {
                        case 'i':
                            table[ln][c] = Integer.parseInt(row[c].trim());
                            break;
                        case 'f':
                            table[ln][c] = Double.parseDouble(row[c].trim());
                            break;
                        case 'c':
                            table[ln][c] = (row[c].trim().replace("'", "")
                                    .replace("\"", "")).charAt(0);
                            break;
                    }
                }
            }

            return table;
        } catch (FileNotFoundException e) {
            System.out.println("Could not find file: '"+path+"'");
        } catch (IOException e) {
            System.out.println("IO Exception!");
        } finally {
            try {
                if (br != null)
                    br.close();
            } catch (IOException e) {
                System.out.println("Couldn't close buffered reader!");
            }
        }
        return new Object[1][1];
    }

    /**
     * A helper function that casts and returns a CSV file consisting
     * entirely of ints on all rows and columns. Useful for loading the
     * atoms table.
     *  @param path The string path to the CSV file.
     *  @return A two dimensional ints array.
     */
    public static int[][] readIntsFile(String path) {
        Object[][] table = readFile(path);
        int rows = table.length;
        int cols = table[0].length;

        int[][] newtable = new int[rows][cols];

        for (int i=0; i<rows; i++)
            for (int j=0; j<cols; j++)
                newtable[i][j] = ((Integer)table[i][j]).intValue();

        return newtable;
    }
}
