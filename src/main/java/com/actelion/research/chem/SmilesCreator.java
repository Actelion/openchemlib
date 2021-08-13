/***************************************************************************
   SmilesCreator.java           (June 05, 1998 / Man-Ling Lee)
Last change:   Jul. 19, 1998

-->Generate a non-unique Smiles string of a given molecule.
-->Stereo features not yet implemented.

Atom description in SMILES:
   [(isotop)(symbol)(implicit H's)(charge)], eg. [15NH2+];
   Implicit H's are needed, if symbol is in a bracket.
   If only the label (of an organic atom) is given, then brackets are not 
      required.
****************************************************************************/

package com.actelion.research.chem;

/**
 * Use the IsomericSmilesCreator instead!!!
 */
@Deprecated
public class SmilesCreator
{
   private ExtendedMolecule mMol;
   protected String mSmiles;

   private int mVisitedMolAtoms;
   private int mMolAtomIsSmiAtom[];     // To get the correspond sAtom.
   private boolean mVisitedMolAtom[];   // The atom order here is same as
   private boolean mVisitedMolBond[];   // the one in the Molecule mol.

   private String mSmiAtomStr[];  // May have another atom order as above.
   private int mRingClosures;     // Number of Ring closures found.
   private int mDisconnections;

   /** 
   public String generateSmiles (Molecule inMol) 
   Effect on instance variable:
       -> increment this.mVisitedMolAtoms
       -> increment this.mDisconnections, 
          if inMol consists of several isolated molecules (or atoms).
       -> Build mSmiles.
   * This class is deprecated. Use the IsomericSmilesCreator instead!!!
   */
   @Deprecated
   public String generateSmiles (ExtendedMolecule inMol)
   {
      int atoms, bonds;
      int i;
      boolean visitedAllAtoms = false;

      // Set variables
      mMol = inMol;

	  mMol.ensureHelperArrays(Molecule.cHelperParities);

      atoms = mMol.getAtoms();
      bonds = mMol.getBonds();

      mVisitedMolBond = new boolean[bonds];
      for (i=0; i<bonds; ++i)
         mVisitedMolBond[i] = false;
      mVisitedMolAtom = new boolean[atoms];
      mMolAtomIsSmiAtom = new int[atoms];
      for (i=0; i<atoms; ++i)
      {
         mVisitedMolAtom[i] = false;
         mMolAtomIsSmiAtom[i] = -1;
      }
      mSmiAtomStr = new String[3*atoms];

      mVisitedMolAtoms = 0;
      mRingClosures = 0;
      mDisconnections = 0;

      // Build the mSmiAtomStr[] for every atoms
      while (visitedAllAtoms == false)
      {
         // Are there any atoms or bonds which are not visited? 
         for (i=0; i<atoms; ++i)
         {
            if (mVisitedMolAtom[i] == false)
            {
               if (mDisconnections > 0)
                  mSmiAtomStr[mVisitedMolAtoms++] = ".";
               visitMolAtom (i, -1);
               ++mDisconnections;
               break;   
            }
         }
         if (i == atoms)
            visitedAllAtoms = true;
      }                                                                                          

      // Build the mSmiles
      mSmiles = "";
      for (i=0; i<mVisitedMolAtoms; ++i)
         mSmiles += mSmiAtomStr[i];
      return(mSmiles);
   }



   /** 
   private void visitMolAtom (int molAtom, int molBond)
   Effect on instance variable:
       -> increment this.mVisitedMolAtoms
       -> increment this.mRingClosures, if rings exist
       -> set this.mSmiAtomStr[mVisitedMolAtoms]
   */
   private void visitMolAtom (int molAtom, int molBond)
   {
      boolean addBracket = true;
      int atomCharge, atomIsotope, atomicNo;
      int branchesToVisit = 0;
      int connAtom, connAtoms, connBond;
      int currentSmiAtom;
      int i, implicitHs;
      int skippedConnAtoms = 0;
      String atomLabel;
      
      currentSmiAtom = mVisitedMolAtoms;
      mMolAtomIsSmiAtom[molAtom] = currentSmiAtom;

      // Get atom properties
      atomicNo = mMol.getAtomicNo(molAtom);
      atomLabel = mMol.getAtomLabel(molAtom);
      atomCharge = mMol.getAtomCharge(molAtom);
      atomIsotope = mMol.getAtomMass(molAtom);
      connAtoms = mMol.getConnAtoms(molAtom);


      // Build mSmiAtomStr[currentSmiAtom]; 
      // The order of the elements are defined by Daylight
      if (atomCharge == 0 && atomIsotope == 0 && isOrganic(atomicNo))
         addBracket = false;
      
      mSmiAtomStr[currentSmiAtom] = "";   // Clear mSmiAtomStr
      if (molBond != -1)
      {
         switch (mMol.getBondOrder(molBond))
         {
            case 0:   mSmiAtomStr[currentSmiAtom] += "~"; break; //Query
            case 2:   mSmiAtomStr[currentSmiAtom] += "="; break;
            case 3:   mSmiAtomStr[currentSmiAtom] += "#"; break;
         }
      }

      if (addBracket == true)
         mSmiAtomStr[currentSmiAtom] += "[";

      if (atomIsotope != 0)
         mSmiAtomStr[currentSmiAtom] += java.lang.Integer.toString(atomIsotope);
      mSmiAtomStr[currentSmiAtom] += atomLabel;
      if (addBracket == true)
      {
         if (0 < (implicitHs = mMol.getImplicitHydrogens(molAtom)))
         {
            mSmiAtomStr[currentSmiAtom] += "H";
            if (1 < implicitHs)
               mSmiAtomStr[currentSmiAtom] += implicitHs;
         }
      }
      if (atomCharge != 0)
      {
         if (atomCharge > 0)
            mSmiAtomStr[currentSmiAtom] += "+";
         else
            mSmiAtomStr[currentSmiAtom] += "-";
         if (Math.abs(atomCharge) > 1)
            mSmiAtomStr[currentSmiAtom] += java.lang.Integer.toString(Math.abs(atomCharge));
      }

      if (addBracket == true)
         mSmiAtomStr[currentSmiAtom] += "]";

      // When ready with current atom, then
      if (molBond != -1)
         mVisitedMolBond[molBond] = true;
      mVisitedMolAtom[molAtom] = true;
      ++mVisitedMolAtoms;


      // Count the branches which has to be visited
      for (i=0; i<connAtoms; ++i)
//       if (mVisitedMolAtom[mMol.getConnAtom(molAtom,i)] == false)		// bug fixed TLS 22.6.2011
         if (mVisitedMolBond[mMol.getConnBond(molAtom,i)] == false)
            ++branchesToVisit;

      // Go to next neigbors
      for (i=0; i<connAtoms; ++i)
      {
         connAtom = mMol.getConnAtom (molAtom, i);
         connBond = mMol.getConnBond (molAtom, i);
         if (mVisitedMolBond[connBond] == true)
         {
            ++skippedConnAtoms; 
            continue;
         }
         if (mVisitedMolAtom[connAtom] == true)   // We have a ring to close.
         {
            ++mRingClosures;
	    mVisitedMolBond[connBond] = true;
            switch (mMol.getBondOrder(connBond))
            {
               case 0:
                  mSmiAtomStr[mMolAtomIsSmiAtom[connAtom]] += "~";
                  mSmiAtomStr[currentSmiAtom] += "~";
                  break;
              case 2:
                  mSmiAtomStr[mMolAtomIsSmiAtom[connAtom]] += "=";
                  mSmiAtomStr[currentSmiAtom] += "=";
                  break;
              case 3:
                  mSmiAtomStr[mMolAtomIsSmiAtom[connAtom]] += "#";
                  mSmiAtomStr[currentSmiAtom] += "3";
                  break;
            }
            if (mRingClosures > 9)
            {
               mSmiAtomStr[mMolAtomIsSmiAtom[connAtom]] += "%";
               mSmiAtomStr[currentSmiAtom] += "%";
            }
            mSmiAtomStr[mMolAtomIsSmiAtom[connAtom]] += java.lang.Integer.toString(mRingClosures);
            mSmiAtomStr[currentSmiAtom] += java.lang.Integer.toString(mRingClosures);
            continue;
         }

		if (i-skippedConnAtoms < branchesToVisit-1)   // new branch
            mSmiAtomStr[mVisitedMolAtoms++] = "(";

         visitMolAtom (connAtom, connBond);

         // end of a branch
         // In meantime all the atoms of the branch starting with 
         // ith atom neighbor of atom currentSmiAtom have been visited.
         // So mVisitedMolAtoms is changed.
         if (i-skippedConnAtoms < branchesToVisit-1)
            mSmiAtomStr[mVisitedMolAtoms++] = ")";
      }
   }


   private boolean isOrganic (int atomicNo)
   {
      switch (atomicNo)
      {
         case  5:                     // B
         case  6:                     // C
         case  7:                     // N
         case  8:                     // O
         case  9:                     // F
         case 15:                     // P
         case 16:                     // S
         case 17:                     // Cl
         case 35:                     // Br
         case 53: return (true);      // I
         default: return (false);     // non-organic atoms in Daylight
      }
   }
}
