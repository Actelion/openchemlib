package com.actelion.research.chem.name;

import com.actelion.research.chem.StereoMolecule;

/**
 * Created by thomas on 7/13/17.
 */
public class StructureNameResolver {
	private static IStructureNameResolver sResolver;

	public static IStructureNameResolver getInstance() {
		return sResolver;
		}

	public static void setInstance(IStructureNameResolver resolver) {
		sResolver = resolver;
	}

	/**
	 * If a IStructureNameResolver instance was instantiated and given to this class, then that is
	 * asked to try to resolve the given chemical name, i.e. to create the respective StereoMolecule
	 * that is represented by the name.
	 * @param name
	 */
	public static StereoMolecule resolve(String name) {
		return sResolver == null ? null : sResolver.resolveName(name);
		}
	}
