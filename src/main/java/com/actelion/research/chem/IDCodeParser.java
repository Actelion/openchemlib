package com.actelion.research.chem;

import com.actelion.research.chem.coords.CoordinateInventor;

public class IDCodeParser extends IDCodeParserWithoutCoordinateInvention {
	private boolean mEnsure2DCoordinates;
	private int mCoordinateMode = CoordinateInventor.MODE_DEFAULT;

	/**
	 * This default constructor creates molecules guaranteed to have 2D-atom-coordinates.
	 * If 2D-coordinates are not supplied with the idcode, or if supplied coordinates are 3D,
	 * then new 2D-coordinates are created on the fly.
	 */
	public IDCodeParser(){
		this(true);
	}

	/**
	 * This default constructor creates molecules guaranteed to have 2D-atom-coordinates.
	 * If 2D-coordinates are not supplied with the idcode, or if supplied coordinates are 3D,
	 * then new 2D-coordinates are created on the fly.
	 * @param coordinateMode mode used for CoordinateInventor
	 */
	public IDCodeParser(int coordinateMode){
		this(true);
		mCoordinateMode = coordinateMode;
	}

	/**
	 * @param ensure2DCoordinates If TRUE and no coordinates are passed with the idcode, then
	 * the parser generates atom coordinates of any molecule and assigns up/down bonds reflecting
	 * given atom parities. Generating coordinates is potentially error-prone, such that providing
	 * original coordinates, where available, should be the preferred option.
	 * <br><b>WARNING:</b> If FALSE: In this case stereo parities are taken directly from the idcode,
	 * missing explicitly 'unknown' parities, because they are not part of the idcode.
	 * Without atom coordinates up/down bonds cannot be assigned. If further processing relies
	 * on up/down bond stereo information or needs to distinguish parities 'none' from 'unknown',
	 * (e.g. idcode creation, checking for stereo centers, calculating the skeletonSpheres descriptor),
	 * or if you are not exactly sure, what to do, then use the constructor IDCodeParser(true).
	 * If you supply encoded 3D-coordinates, then use IDCodeParser(false).
	 */
	public IDCodeParser(boolean ensure2DCoordinates) {
		super();
		mEnsure2DCoordinates = ensure2DCoordinates;
		}



	@Override
	protected boolean ensure2DCoordinates() {
		return mEnsure2DCoordinates;
		}

	@Override
	protected void inventCoordinates(StereoMolecule mol) {
		CoordinateInventor inventor = new CoordinateInventor(mCoordinateMode);
		inventor.setRandomSeed(0x1234567890L);  // create reproducible coordinates
		inventor.invent(mol);
		}
	}
