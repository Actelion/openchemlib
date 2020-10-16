package com.actelion.research.chem.docking;

public class DockingFailedException extends Exception { 
    public DockingFailedException(String errorMessage) {
        super(errorMessage);
    }
}