package com.actelion.research.chem.name;

import com.actelion.research.chem.StereoMolecule;

/**
 * Created by thomas on 7/13/17.
 */
public interface IStructureNameResolver {
	/**
	 * Local and typically quick name resolution
 	 */
	public StereoMolecule resolveLocal(String name);

	/**
	 * Typically remote server based name resolution that requires a network round trip
	 */
	public StereoMolecule resolveRemote(String name);

	/**
	 * Remote server based name resolution of multiple names packed into one network round trip
	 */
	public String[] resolveRemote(String[] nameList);
}
