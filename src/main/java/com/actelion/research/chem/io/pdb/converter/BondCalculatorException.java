package com.actelion.research.chem.io.pdb.converter;
/*
 * This file is part of cif2molecule.
 *
 * cif2molecule is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * @author Antanas Vaitkus
 */

/**
 * An exception class to be used with the BondCalculator class.
 */
class BondCalculatorException extends Exception {
    BondCalculatorException(String message) {
        super(message);
    }
}
