/**
 * TODO: Add license header (switch to MIT)
 */
package de.unijena.cheminf.deglycosylation;

/**
 * TODO:
 * - Is there a way to combine overlapping matches of linear patterns? Would that be better? Maybe do not remove linear
 * sugar matches if they are too small after removal of overlapping atoms, e.g. only a CH3 group could get removed (see tests)
 * - Check the exceptions the private methods throw, are they passed on by the public methods?
 * - Add MIT license header (and change license on GitHub)
 * - make a console application to deploy as jar; see google doc for input and output
 * - add more linear sugars, e.g. aldopentose, different deoxypentoses/deoxyhexoses, heptoses, octoses
 * - see to dos in the code
 * - check occurrence of 'questionable' linear sugars
 * - return a list of the removed sugars; that will be complicated if the given molecule should be cloned!
 * - maybe make some more methods public that might be of interest
 *
 * To discuss:
 * - see 'to do / to discuss' points in the code
 */

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerComparator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utility class to remove sugar moieties from molecular structures, primarily natural products.
 *
 * @author Jonas Schaub, Maria Sorokina
 */
public class SugarRemovalUtility {
    //<editor-fold desc="Enum StructuresToKeepMode">
    /**
     * Enum that contains options for how to judge whether a remaining (unconnected) substructure after sugar removal is
     * worth keeping or should be discarded because it is too small/light etc.
     * <br>The set option also plays a crucial role in judging whether a sugar moiety is terminal or not.
     * <br>Also, the set threshold for the molecular weight / heavy atom count etc. interrelates with this option.
     */
    public static enum StructuresToKeepMode {
        /**
         * Specifies that all structures are worth keeping.
         */
        ALL (0),

        /**
         * Specifies that whether a structure is worth keeping will be judged by its heavy atom count. The default
         * threshold to keep a structure is set to 5 heavy atoms.
         */
        HEAVY_ATOM_COUNT (5),

        /**
         * Specifies that whether a structure is worth keeping will be judged by its molecular weight. The default
         * threshold to keep a structure is set to 60 Da (= 5 carbon atoms).
         */
        MOLECULAR_WEIGHT (60);

        /**
         * Default threshold to keep a structure for the respective option.
         */
        private final int defaultThreshold;

        /**
         * Constructor.
         *
         * @param aDefaultValue the default threshold to keep a structure for the respective option; no parameter checks
         *                      are implemented
         */
        StructuresToKeepMode(int aDefaultValue) {
            this.defaultThreshold = aDefaultValue;
        }

        /**
         * Returns the default threshold to keep a structure for this option
         *
         * @return the default threshold
         */
        public int getDefaultThreshold() {
            return this.defaultThreshold;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * Property key to indicate that the structure contains (or contained before removal) circular sugar moieties.
     */
    public static final String CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY = "CONTAINS_CIRCULAR_SUGAR";

    /**
     * Property key to indicate that the structure contains (or contained before removal) linear sugar moieties.
     */
    public static final String CONTAINS_LINEAR_SUGAR_PROPERTY_KEY = "CONTAINS_LINEAR_SUGAR";

    /**
     * Property key to indicate that the structure contains (or contained before removal) sugar moieties.
     */
    public static final String CONTAINS_SUGAR_PROPERTY_KEY = "CONTAINS_SUGAR";

    /**
     * Property key for index that is added to any IAtom object in a given IAtomContainer object for internal unique
     * identification of the respective IAtom object.
     */
    public static final String INDEX_PROPERTY_KEY = "SUGAR_REMOVAL_UTILITY_INDEX";

    /**
     * Linear sugar structures represented as SMILES codes. An input molecule is scanned for these substructures for
     * the detection of linear sugars.
     */
    public static final String[] LINEAR_SUGARS_SMILES = {
            "C(C(C(C(C(C=O)O)O)O)O)O", //aldohexose
            "C(C(CC(C(CO)O)O)O)(O)=O", //3-deoxyhexonic acid
            "C(C(C(CC(=O)O)O)O)O", //2-deoxypentonic acid
            "C(C(C(C(C(CO)O)O)O)=O)O", //2-ketohexose
            "C(C(C(C(C(CO)O)O)O)O)O", //hexitol
            "C(C(C(C(CC=O)O)O)O)O", //2-deoxyhexose
            "CCCCC(O)C(=O)O", //2-hydroxyhexanoic acid TODO/discuss: Is this a sugar?
            "CC(=O)CC(=O)CCC(=O)O", //4,6-dioxoheptanoic acid TODO/discuss: Is this a sugar?
            "O=C(O)CC(O)CC(=O)O", //3-hydroxypentanedioic acid TODO/discuss: Is this a sugar?
            "O=C(O)C(=O)C(=O)C(O)C(O)CO", //hexo-2,3-diulosonic acid
            "O=C(O)CCC(O)C(=O)O", //2-hydroxypentanedioic acid TODO/discuss: Is this a sugar?
            "O=C(CO)C(O)C(O)CO" //ketopentose
    };

    /**
     * Circular sugar structures represented as SMILES codes. The isolated rings of an input molecule are matched with
     * these structures for the detection of circular sugars. The structures listed here only represent the circular
     * part of sugar rings (i.e. one oxygen atom and multiple carbon atoms). Common exocyclic structures like
     * hydroxy groups are not part of the patterns and detected in another step.
     */
    public static final String [] RING_SUGARS_SMILES = {
            "C1CCOC1", //tetrahydrofuran to match all 5-membered sugar rings
            "C1CCOCC1", //tetrahydropyran to match all 6-membered sugar rings
            "C1CCCOCC1" //oxepane to match all 7-membered sugar rings
    };

    /**
     * Default setting for whether glycosidic bonds between sugar rings should be detected to determine whether a
     * candidate sugar structure should be removed (default: false).
     */
    public static final boolean DETECT_GLYCOSIDIC_BOND_DEFAULT = false;

    /**
     * Default setting for whether only terminal sugar moieties should be removed, i.e. those that result in a still
     * fully-connected structure (default: true).
     */
    public static final boolean REMOVE_ONLY_TERMINAL_DEFAULT = true;

    /**
     * Default setting for how to judge whether a remaining (unconnected) substructure after sugar removal is
     * worth keeping or should be discarded because it is too small/light etc (default: judge by heavy atom count).
     * The minimum value to reach for the respective characteristic to judge by is set in an additional option and all
     * enum constants have their own default values. See the StructureToKeepMode enum.
     */
    public static final StructuresToKeepMode STRUCTURES_TO_KEEP_MODE_DEFAULT = StructuresToKeepMode.HEAVY_ATOM_COUNT;

    /**
     * Default setting for whether the number of attached, exocyclic, single-bonded oxygen atoms should be evaluated to
     * determine whether a circular candidate sugar structure should be removed (default: true).
     */
    public static final boolean INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT = true;

    /**
     * Default setting for the minimum ratio of attached, exocyclic, single-bonded oxygen atoms to the number of atoms
     * in the candidate circular sugar structure to reach in order to be classified as a sugar moiety
     * if the number of exocyclic oxygen atoms should evaluated (default: 0.5 so at a minimum 3 connected, exocyclic
     * oxygen atoms for a six-membered ring).
     */
    public static final double ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT = 0.5;

    /**
     * Default setting for whether linear sugar structures that are part of a larger ring should be removed (default:
     * false).
     */
    public static final boolean REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT = false;

    /**
     * Default setting for whether to add a property to given atom containers to indicate that the structure contains
     * (or contained before removal) sugar moieties. See property keys in the public constants of this class.
     */
    public static final boolean SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT = true;
    //</editor-fold>
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(SugarRemovalUtility.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * Universal isomorphism tester for the detection of circular sugars.
     */
    private final UniversalIsomorphismTester univIsomorphismTester;
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    /**
     * Linear sugar structures parsed into atom containers. Not used for detection but parsed into patterns.
     */
    private List<IAtomContainer> linearSugars;

    /**
     * Circular sugar structures parsed into atom containers. Used for detection via universal isomorphism tester.
     */
    private List<IAtomContainer> ringSugars;

    /**
     * Patterns of linear sugar structures to detect linear sugar moieties in the given molecules.
     */
    private List<DfPattern> linearSugarPatterns;

    /**
     * Detect glycosidic bond setting.
     */
    private boolean detectGlycosidicBond;

    /**
     * Remove only terminal sugar moieties setting.
     */
    private boolean removeOnlyTerminal;

    /**
     * Structure to keep mode setting (see enum at the top).
     */
    private StructuresToKeepMode structuresToKeepMode;

    /**
     * Threshold for the characteristic of an unconnected fragment set in the structure to keep mode to judge whether to
     * keep it or discard it.
     */
    private int structureToKeepModeThreshold;

    /**
     * Include number/ratio of connected, exocyclic oxygen atoms setting.
     */
    private boolean includeNrOfAttachedOxygens;

    /**
     * Minimum ratio of connected, exocyclic oxygen atoms to the number of atoms in the candidate sugar ring.
     */
    private double attachedOxygensToAtomsInRingRatioThreshold;

    /**
     * Remove linear sugars in circular structures setting.
     */
    private boolean removeLinearSugarsInRing;

    /**
     * Add a property to sugar-containing, given atom containers setting.
     */
    private boolean setPropertyOfSugarContainingMolecules;
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * Sole constructor of this class. The circular and linear sugar structures are parsed into atom containers or pattern
     * and all settings are set to their default values (see public static constants or enquire via get/is methods).
     */
    public SugarRemovalUtility() {
        this.univIsomorphismTester = new UniversalIsomorphismTester();
        this.linearSugars = new ArrayList<>(SugarRemovalUtility.LINEAR_SUGARS_SMILES.length);
        this.ringSugars = new ArrayList<>(SugarRemovalUtility.RING_SUGARS_SMILES.length);
        this.linearSugarPatterns = new ArrayList<>(SugarRemovalUtility.LINEAR_SUGARS_SMILES.length);
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //adding linear sugars to list
        for (String tmpSmiles : SugarRemovalUtility.LINEAR_SUGARS_SMILES) {
            try {
                this.linearSugars.add(tmpSmilesParser.parseSmiles(tmpSmiles));
            } catch (InvalidSmilesException anInvalidSmilesException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, anInvalidSmilesException.toString(), anInvalidSmilesException);
            }
        }
        //adding ring sugars to list
        for (String tmpSmiles : SugarRemovalUtility.RING_SUGARS_SMILES) {
            try {
                this.ringSugars.add(tmpSmilesParser.parseSmiles(tmpSmiles));
            } catch (InvalidSmilesException anInvalidSmilesException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, anInvalidSmilesException.toString(), anInvalidSmilesException);
            }
        }
        //parsing linear sugars into patterns
        for(IAtomContainer tmpSugarAC : this.linearSugars){
            this.linearSugarPatterns.add(DfPattern.findSubstructure(tmpSugarAC));
        }
        this.detectGlycosidicBond = SugarRemovalUtility.DETECT_GLYCOSIDIC_BOND_DEFAULT;
        this.removeOnlyTerminal = SugarRemovalUtility.REMOVE_ONLY_TERMINAL_DEFAULT;
        this.structuresToKeepMode = SugarRemovalUtility.STRUCTURES_TO_KEEP_MODE_DEFAULT;
        this.structureToKeepModeThreshold = this.structuresToKeepMode.defaultThreshold;
        this.includeNrOfAttachedOxygens = SugarRemovalUtility.INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT;
        this.attachedOxygensToAtomsInRingRatioThreshold =
                SugarRemovalUtility.ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT;
        this.removeLinearSugarsInRing = SugarRemovalUtility.REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT;
        this.setPropertyOfSugarContainingMolecules = SugarRemovalUtility.SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get/is">
    /**
     * Get a list of (unique) SMILES strings representing the linear sugar structures an input molecule is scanned for
     * by the methods for sugar detection or removal. The default structures can also be retrieved from the respective
     * public constant of this class. But additional structures added externally are only returned by this method.
     * <br>If a structure cannot be parsed into a SMILES string, it is excluded from the list.
     *
     * @return a list of SMILES codes
     */
    public List<String> getLinearSugars() {
        List<String> tmpSmilesList = new ArrayList<>(this.linearSugars.size());
        SmilesGenerator tmpSmilesGen = new SmilesGenerator(SmiFlavor.Unique);
        for (IAtomContainer tmpLinearSugar : this.linearSugars) {
            String tmpSmiles = null;
            try {
                tmpSmiles = tmpSmilesGen.create(tmpLinearSugar);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
            }
            if (!Objects.isNull(tmpSmiles)) {
                tmpSmilesList.add(tmpSmiles);
            }
        }
        return tmpSmilesList;
    }

    /**
     * Get a list of (unique) SMILES strings representing the circular sugar structures an input molecule is scanned for
     * by the methods for sugar detection or removal. The default structures can also be retrieved from the respective
     * public constant of this class. But additional structures added externally are only returned by this method.
     * <br>If a structure cannot be parsed into a SMILES string, it is excluded from the list.
     *
     * @return a list of SMILES codes
     */
    public List<String> getCircularSugars() {
        List<String> tmpSmilesList = new ArrayList<>(this.ringSugars.size());
        SmilesGenerator tmpSmilesGen = new SmilesGenerator(SmiFlavor.Unique);
        for (IAtomContainer tmpRingSugar : this.ringSugars) {
            String tmpSmiles = null;
            try {
                tmpSmiles = tmpSmilesGen.create(tmpRingSugar);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
            }
            if (!Objects.isNull(tmpSmiles)) {
                tmpSmilesList.add(tmpSmiles);
            }
        }
        return tmpSmilesList;
    }

    /**
     * Specifies whether the glycosidic bond of a circular candidate sugar structure is detected and its presence taken
     * into account to decide on whether the sugar is to be removed.
     *
     * @return true if only circular sugar moieties connected via a glycosidic bond are removed according to the current
     * settings
     */
    public boolean isGlycosidicBondDetected() {
        return this.detectGlycosidicBond;
    }

    /**
     * Specifies whether only terminal sugar moieties (i.e. those that result in a still fully-connected structure) are
     * removed.
     *
     * @return true if only terminal sugar moieties are removed according to the current settings
     */
    public boolean areOnlyTerminalSugarsRemoved() {
        return this.removeOnlyTerminal;
    }

    /**
     * Returns the current setting on how to judge whether an unconnected substructure resulting from the sugar removal
     * should be discarded or kept. This can e.g. be judged by its heavy atom count or its molecular weight or it can be
     * set that all structures are to be kept. If too small / too light structures are discarded, an additional threshold
     * is specified that the structures have to fulfill in order to be kept (i.e. to be judged 'big/heavy enough').
     *
     * @return a StructuresToKeepMode enum object representing the current setting
     */
    public StructuresToKeepMode getStructuresToKeepMode() {
        return this.structuresToKeepMode;
    }

    /**
     * Returns the current threshold of e.g. molecular weight or heavy atom count (depending on the currently set
     * structure to keep mode) an unconnected substructure resulting from the sugar removal has to reach in order to be
     * kept and not discarded.
     *
     * @return an integer specifying the currently set threshold
     */
    public int getStructureToKeepModeThreshold() {
        return this.structureToKeepModeThreshold;
    }

    /**
     * Specifies whether the number of attached, exocyclic, single-bonded oxygen atoms is currently evaluated to determine
     * whether a circular candidate sugar structure is to be removed. If this option is set, the circular sugar candidates
     * have to reach an additionally specified minimum ratio of said oxygen atoms to the number of atoms in the respective
     * ring in order to be seen as a sugar ring and therefore removed.
     *
     * @return true if the ratio of the number of attached, exocyclic, single-bonded oxygen atoms to the number of atoms
     * in the candidate sugar ring is included in the decision making process according to the current settings
     */
    public boolean isNrOfAttachedOxygensIncluded() {
        return this.includeNrOfAttachedOxygens;
    }

    /**
     * Returns the currently set minimum ratio of attached, exocyclic single-bonded oxygen atoms to the number of atoms
     * in the candidate circular sugar structure to be classified as a sugar moiety and therefore removed
     * (if the option to include this characteristic in the decision making process is enabled).
     *
     * @return the minimum ratio of attached oxygen atoms to the number of atoms in the sugar ring
     */
    public double getAttachedOxygensToAtomsInRingRatioThreshold() {
        return this.attachedOxygensToAtomsInRingRatioThreshold;
    }

    /**
     * Specifies whether linear sugar structures that are part of a larger ring are removed according to the current
     * settings.
     *
     * @return true if linear sugars in larger rings are removed with the current settings
     */
    public boolean areLinearSugarsInRingsRemoved() {
        return this.removeLinearSugarsInRing;
    }

    /**
     * Specifies whether a respective property is added to given atom containers that contain (or contained before
     * removal) sugar moieties. See property keys in the public constants of this class.
     *
     * @return true if properties are added to the given atom containers
     */
    public boolean isPropertyOfSugarContainingMoleculesSet() {
        return this.setPropertyOfSugarContainingMolecules;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set/add/clear">
    /**
     * Allows to add an additional sugar ring to the list of circular sugar structures an input molecule is scanned for
     * by the methods for sugar detection and removal. The given structure must not be isomorph to the already present
     * ones and it must contain only exactly one isolated ring (because only the isolated rings of an input structure are
     * matched with the circular sugar pattern).
     *
     * @param aCircularSugar an atom container representing only one isolated sugar ring
     * @throws NullPointerException if given atom container is 'null'
     * @throws IllegalArgumentException if the given atom container is empty or does represent a molecule that contains
     * no isolated ring, more than one isolated ring, consists of more structures than one isolated ring or is isomorph
     * to a circular sugar structure already present
     */
    public void addCircularSugar(IAtomContainer aCircularSugar) throws NullPointerException, IllegalArgumentException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aCircularSugar, "Given atom container is 'null'");
        if (aCircularSugar.isEmpty()) {
            throw new IllegalArgumentException("Given atom container is empty.");
        }
        int[][] tmpAdjList = GraphUtil.toAdjList(aCircularSugar);
        RingSearch tmpRingSearch = new RingSearch(aCircularSugar, tmpAdjList);
        List<IAtomContainer> tmpIsolatedRingFragments = tmpRingSearch.isolatedRingFragments();
        int tmpSize = tmpIsolatedRingFragments.size();
        if (!(tmpSize == 1)) {
            throw new IllegalArgumentException("Given molecule is either not circular or contains too many rings.");
        }
        UniversalIsomorphismTester tmpUnivIsomorphTester = new UniversalIsomorphismTester();
        boolean tmpIsolatedRingMatchesEntireInputStructure = false;
        IAtomContainer tmpIsolatedRing = tmpIsolatedRingFragments.get(0);
        try {
            tmpIsolatedRingMatchesEntireInputStructure = tmpUnivIsomorphTester.isIsomorph(aCircularSugar, tmpIsolatedRing);
        } catch (CDKException aCDKException) {
            SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
        }
        if (!tmpIsolatedRingMatchesEntireInputStructure) {
            throw new IllegalArgumentException("The given structure does not only consist of one isolated ring.");
        }
        for (IAtomContainer tmpSugar : this.ringSugars) {
            boolean tmpIsIsomorph = false;
            try {
                tmpIsIsomorph = tmpUnivIsomorphTester.isIsomorph(tmpSugar, aCircularSugar);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
                throw new IllegalArgumentException("Could not determine isomorphism with already present sugar structures.");
            }
            if (tmpIsIsomorph) {
                throw new IllegalArgumentException("Given sugar pattern is already present.");
            }
        }
        //</editor-fold>
        this.ringSugars.add(aCircularSugar);
    }

    /**
     * Allows to add an additional sugar ring (represented as a SMILES string) to the list of circular sugar structures
     * an input molecule is scanned for by the methods for sugar detection and removal. The given structure must not be
     * isomorph to the already present ones and it must contain only exactly one isolated ring (because only the isolated
     * rings of an input structure are matched with the circular sugar pattern).
     *
     * @param aSmilesCode a SMILES code representation of a molecule consisting of only one isolated sugar ring
     * @throws NullPointerException if given string is 'null'
     * @throws IllegalArgumentException if the given SMILES string is empty or does represent a molecule that contains
     * no isolated ring, more than one isolated ring, consists of more structures than one isolated ring or is isomorph
     * to a circular sugar structure already present
     * @throws CDKException if the given SMILES string cannot be parsed into a molecular structure
     */
    public void addCircularSugar(String aSmilesCode) throws NullPointerException, CDKException, IllegalArgumentException {
        Objects.requireNonNull(aSmilesCode, "Given SMILES code is 'null'");
        if (aSmilesCode.isEmpty()) {
            throw new IllegalArgumentException("Given SMILES code is empty");
        }
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpRingSugar = tmpSmiPar.parseSmiles(aSmilesCode);
        //see exceptions this method can throw!
        this.addCircularSugar(tmpRingSugar);
    }

    /**
     * Allows to add an additional linear sugar to the list of linear sugar structures an input molecule is scanned for
     * by the methods for sugar detection and removal. The given structure must not be isomorph to the already present
     * ones (no further requirements, so in fact, any kind of structure could be added here, e.g. to detect/remove amino
     * acids also). Note: If the given structure contains circles, to remove its matches entirely in the removal methods,
     * the option to remove linear sugars in rings needs to be disabled. Otherwise, all circular substructures of the
     * 'linear sugars' will not be removed.
     *
     *
     * @param aLinearSugar an atom container representing a molecular structure to search for
     * @throws NullPointerException if given atom container is 'null'
     * @throws IllegalArgumentException if the given atom container is empty or is isomorph to a linear sugar structure
     * already present
     */
    public void addLinearSugar(IAtomContainer aLinearSugar) throws NullPointerException, IllegalArgumentException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aLinearSugar, "Given atom container is 'null'");
        if (aLinearSugar.isEmpty()) {
            throw new IllegalArgumentException("Given atom container is empty.");
        }
        //note: no check for linearity here to allow adding of structures that contain rings, e.g. amino acids
        UniversalIsomorphismTester tmpUnivIsomorphTester = new UniversalIsomorphismTester();
        for (IAtomContainer tmpSugar : this.linearSugars) {
            boolean tmpIsIsomorph = false;
            try {
                tmpIsIsomorph = tmpUnivIsomorphTester.isIsomorph(tmpSugar, aLinearSugar);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
                throw new IllegalArgumentException("Could not determine isomorphism with already present sugar structures.");
            }
            if (tmpIsIsomorph) {
                throw new IllegalArgumentException("Given sugar pattern is already present.");
            }
        }
        //</editor-fold>
        this.linearSugars.add(aLinearSugar);
        this.linearSugarPatterns.add(DfPattern.findSubstructure(aLinearSugar));
    }

    /**
     * Allows to add an additional linear sugar (represented as a SMILES string) to the list of linear sugar structures
     * an input molecule is scanned for by the methods for sugar detection and removal. The given structure must not be
     * isomorph to the already present ones (no further requirements, so in fact, any kind of structure could be added
     * here, e.g. to detect/remove amino acids also). Note: If the given structure contains circles, to remove its
     * matches entirely in the removal methods, the option to remove linear sugars in rings needs to be disabled.
     * Otherwise, all circular substructures of the 'linear sugars' will not be removed.
     *
     * @param aSmilesCode a SMILES code representation of a molecular structure to search for
     * @throws NullPointerException if given string is 'null'
     * @throws IllegalArgumentException if the given SMILES string is empty or does represent a molecule that is isomorph
     * to a linear sugar structure already present
     * @throws CDKException if the given SMILES string cannot be parsed into a molecular structure
     */
    public void addLinearSugar(String aSmilesCode) throws NullPointerException, CDKException, IllegalArgumentException {
        Objects.requireNonNull(aSmilesCode, "Given SMILES code is 'null'");
        if (aSmilesCode.isEmpty()) {
            throw new IllegalArgumentException("Given SMILES code is empty");
        }
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpLinearSugar = tmpSmiPar.parseSmiles(aSmilesCode);
        //see exceptions this method can throw!
        this.addLinearSugar(tmpLinearSugar);
    }

    /**
     * Internally clears all the circular sugar structures an input molecule is scanned for by the methods for sugar
     * detection and removal.
     */
    public void clearCircularSugars() {
        this.ringSugars.clear();
    }

    /**
     * Internally clears all the linear sugar structures an input molecule is scanned for by the methods for sugar
     * detection and removal.
     */
    public void clearLinearSugars() {
        this.linearSugars.clear();
        this.linearSugarPatterns.clear();
    }

    /**
     * Sets the option to detect the glycosidic bond of a circular candidate sugar structure and to remove only those
     * having such a bond.
     *
     * @param aBoolean true, if only circular sugar moieties connected via a glycosidic bond should be removed
     */
    public void setDetectGlycosidicBond(boolean aBoolean) {
        this.detectGlycosidicBond = aBoolean;
    }

    /**
     * Sets the option to remove only terminal sugar moieties (i.e. those that result in a still fully-connected
     * structure).
     *
     * @param aBoolean true, if only terminal sugar moieties should be removed
     */
    public void setRemoveOnlyTerminalSugars(boolean aBoolean) {
        this.removeOnlyTerminal = aBoolean;
    }

    /**
     * Sets the characteristic an unconnected substructure resulting from the sugar removal is judged by to determine
     * whether it is worth keeping or should be discarded. This characteristic can e.g. be its heavy atom count or its
     * molecular weight. Alternatively, it can be specified that all structures are to be kept. The available options
     * can be selected from the StructureToKeepMode enum. If too small / too light structures are discarded, an
     * additional threshold is specified that the structures have to fulfill in order to be kept (i.e. to be judged
     * 'big/heavy enough').
     * <br>IMPORTANT NOTE: This threshold is overridden with the default threshold for the newly set option if this
     * method is called!
     * <br>Note: If the structure to keep mode 'All' is combined with the removal of only terminal sugars, even though
     * seemingly terminal circular sugars are not removed if they would leave behind their hydroxy groups as unconnected
     * fragments upon removal.
     *
     * @param aMode the selected mode
     * @throws NullPointerException if given mode is 'null'
     */
    public void setStructuresToKeepMode(StructuresToKeepMode aMode) throws NullPointerException {
        Objects.requireNonNull(aMode, "Given mode is 'null'.");
        this.structuresToKeepMode = aMode;
        this.structureToKeepModeThreshold = this.structuresToKeepMode.getDefaultThreshold();
    }

    /**
     * Sets the threshold of e.g. molecular weight or heavy atom count (depending on the currently set structure to keep
     * mode) an unconnected substructure resulting from the sugar removal has to reach in order to be kept and not
     * discarded.
     *
     * @param aThreshold the new threshold
     * @throws IllegalArgumentException if the structure to keep mode is currently set to keep all structures or the
     * threshold is negative
     */
    public void setStructuresToKeepThreshold(int aThreshold) throws IllegalArgumentException {
        if ((this.structuresToKeepMode == StructuresToKeepMode.ALL)) {
            throw new IllegalArgumentException("The mode is currently set to keep all structures, so a threshold " +
                    "makes no sense.");
        }
        if (aThreshold < 0) {
            throw new IllegalArgumentException("Threshold cannot be negative.");
        }
        this.structureToKeepModeThreshold = aThreshold;
    }

    /**
     * Sets the option to evaluate the number of attached, exocyclic, single-bonded oxygen atoms to determine
     * whether a circular candidate sugar structure is to be removed. If this option is set, the circular sugar candidates
     * have to reach an additionally specified minimum ratio of said oxygen atoms to the number of atoms in the respective
     * ring in order to be seen as a sugar ring and therefore removed.
     *
     * @param aBoolean true, if only sugar rings with a sufficient number of connected oxygen atoms should be removed
     */
    public void setIncludeNrOfAttachedOxygens(boolean aBoolean) {
        this.includeNrOfAttachedOxygens = aBoolean;
    }

    /**
     * Sets the minimum ratio of attached, exocyclic single-bonded oxygen atoms to the number of atoms
     * in the candidate circular sugar structure to be classified as a sugar moiety and therefore removed
     * (if the option to include this characteristic in the decision making process is enabled).
     * <br>E.g.: a ratio of 0.5 means that a six-membered candidate sugar ring needs to have at least 3 attached, exocyclic
     * single-bonded oxygen atoms in order to classify as a sugar.
     * <br>Note: The normally present oxygen atom within a sugar ring is included in the number of ring atoms. So setting
     * the threshold to 1.0 implies that at least one of the carbon atoms in the ring has two attached oxygen atoms.
     * In general, the threshold can be set to values higher than 1.0 but it does not make a lot of sense.
     *
     * @param aDouble the new ratio threshold
     * @throws IllegalArgumentException if the given number is infinite, 'NaN' or smaller than 0 or if the ratio is not
     * evaluated under the current settings
     */
    public void setAttachedOxygensToAtomsInRingRatioThreshold(double aDouble) throws IllegalArgumentException {
        //<editor-fold desc="Checks">
        boolean tmpIsFinite = Double.isFinite(aDouble); //false for NaN and infinity arguments
        boolean tmpIsNegative = (aDouble < 0);
        if(!tmpIsFinite || tmpIsNegative) {
            throw new IllegalArgumentException("Given double is NaN, infinite or negative.");
        }
        if ((!this.includeNrOfAttachedOxygens)) {
            throw new IllegalArgumentException("The number of attached oxygen atoms is currently not included in the " +
                    "decision making process, so a ratio threshold makes no sense.");
        }
        //</editor-fold>
        this.attachedOxygensToAtomsInRingRatioThreshold = aDouble;
    }

    /**
     * Sets the option to remove linear sugar structures even if they are part of a larger ring.
     *
     * @param aBoolean true, if the linear sugar structures that are part of a ring should be removed
     */
    public void setRemoveLinearSugarsInRing(boolean aBoolean) {
        this.removeLinearSugarsInRing = aBoolean;
    }

    /**
     * Sets the option to add a respective property to given atom containers that contain (or contained before
     * removal) sugar moieties. See property keys in the public constants of this class.
     *
     * @param aBoolean true, if properties should be added to the given atom containers
     */
    public void setPropertyOfSugarContainingMolecules(boolean aBoolean) {
        this.setPropertyOfSugarContainingMolecules = aBoolean;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Scans the given molecule for the presence of linear sugars. The result of this method is only influenced by the
     * setting specifying whether linear sugars that are part of ring structures should be removed. It is not influenced
     * by e.g. the setting specifying whether only terminal sugar moieties should be removed.
     * <br>If the respective option is set, a property will be added to the given atom container specifying whether
     * it contains (linear) sugar moieties or not (in addition to the return value of this method).
     *
     * @param aMolecule the atom container to scan for the presence of linear sugar moieties
     * @return true, if the given molecule contains linear sugar moieties
     * @throws NullPointerException if the given atom container is 'null'
     */
    public boolean hasLinearSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        aMolecule = this.setIndices(aMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getLinearSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        return tmpContainsSugar;
    }

    /**
     * Scans the given molecule for the presence of circular sugars. The result of this method is influenced by the
     * settings specifying whether only sugar rings with a glycosidic bond should be removed and whether the number of
     * attached, exocyclic single-bonded oxygen atoms of a candidate sugar ring need to be of sufficient number to classify
     * the ring as a sugar. It is not influenced e.g. the setting specifying whether only terminal sugar moieties should
     * be removed.
     * <br>If the respective option is set, a property will be added to the given atom container specifying whether
     * it contains (circular) sugar moieties or not (in addition to the return value of this method).
     *
     * @param aMolecule the atom container to scan for the presence of circular sugar moieties
     * @return true, if the given molecule contains circular sugar moieties
     * @throws NullPointerException if the given atom container is 'null'
     */
    public boolean hasCircularSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        aMolecule = this.setIndices(aMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        return tmpContainsSugar;
    }

    /**
     * Scans the given molecule for the presence of linear and circular sugars. The result of this method is influenced
     * by some of the settings, e.g. the one specifying whether linear sugars that are part of ring structures should be
     * removed. It is not influenced by e.g. the setting specifying whether only terminal sugar moieties should be removed.
     * <br>If the respective option is set, a property will be added to the given atom container specifying whether
     * it contains (linear/circular) sugar moieties or not (in addition to the return value of this method).
     *
     * @param aMolecule the atom container to scan for the presence of sugar moieties
     * @return true, if the given molecule contains sugar moieties
     * @throws NullPointerException if the given atom container is 'null'
     */
    public boolean hasSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        aMolecule = this.setIndices(aMolecule);
        List<IAtomContainer> tmpCircularSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsCircularSugar = !tmpCircularSugarCandidates.isEmpty();
        List<IAtomContainer> tmpLinearSugarCandidates = this.getLinearSugarCandidates(aMolecule);
        boolean tmpContainsLinearSugar = !tmpLinearSugarCandidates.isEmpty();
        boolean tmpContainsSugar = (tmpContainsCircularSugar || tmpContainsLinearSugar);
        if (this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsCircularSugar);
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsLinearSugar);
        }
        return tmpContainsSugar;
    }

    /**
     * Removes circular sugar moieties from the given atom container. Which substructures are removed depends on the
     * various settings available in this class.
     * <p>
     *     First, indices for a unique identification are added to every IAtom object of the given atom container object
     *     as a property.
     *     <br>Then, the isolated rings of the molecule are matched with the internal circular sugar structure patterns.
     *     All matching rings are candidate sugar structures.
     *     <br>Discarded from the candidates are rings that have an exocyclic double or triple bond. If the respective
     *     option is set, all candidates are discarded that do not have a glycosidic bond. Furthermore, all
     *     candidates are discarded that have too few attached exocyclic oxygen atoms if the respective option is set.
     *     <br>If only terminal sugar moieties are to be removed, the sugar candidates are one-by-one tested for
     *     whether they are terminal or not and removed if they are. The procedure starts anew after iterating over all
     *     candidates and stops if no terminal sugar was removed in one whole iteration. The determination whether
     *     a sugar is terminal or not depends also on the set structure to keep mode because if the removal of a sugar
     *     candidate leads to an unconnected fragment the sugar candidate can still be terminal if the fragment is
     *     too small to keep and therefore discarded after the removal of the sugar. Resulting again in one connected
     *     remaining structure.
     *     <br>If also non-terminal sugars should be removed, simply all detected candidates are removed from the given
     *     molecular structure. If too small/light structures are also removed according to the current settings, the
     *     resulting structure will be cleared of all unconnected fragments that are not big/heavy enough. Still, a
     *     molecule consisting of multiple unconnected structures may result.
     *     <br>If the respective option is set, a property will be added to the given atom container specifying whether
     *     it contains (circular) sugar moieties or not.
     * </p>
     *
     * @param aMolecule the molecule to remove circular sugar moieties from
     * @param aShouldBeCloned true, if the sugar moieties should not be removed from the given atom container but a clone
     *                        of it should be generated and the sugars be removed from that
     * @return if the given atom container should NOT be cloned, this method returns the same given atom container after the
     * sugar removal; the returned molecule may be unconnected if also non-terminal sugars are removed according to
     * the settings and it may be empty if the resulting structure after sugar removal was too small to keep due to the
     * set structure to keep mode and the associated threshold (i.e. the molecule basically was a sugar)
     * @throws NullPointerException if the given atom container is 'null'
     * @throws CloneNotSupportedException if the given atom container should be cloned but does not allow cloning
     * @throws IllegalArgumentException if only terminal sugars should be removed but the given atom container already
     * contains multiple, unconnected structures which makes the determination of terminal and non-terminal structures
     * impossible
     */
    public IAtomContainer removeCircularSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures.");
            }
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        tmpNewMolecule = this.setIndices(tmpNewMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getCircularSugarCandidates(tmpNewMolecule);
        /*note: this means that there are matches of the circular sugar patterns and that they adhere to most of
        the given settings. The exception is that they might not be terminal*/
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        if (tmpContainsSugar) {
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpSugarCandidates);
            tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        }
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    /**
     * Removes linear sugar moieties from the given atom container. Which substructures are removed depends on the
     * various settings available in this class.
     * <p>
     *     First, indices for a unique identification are added to every IAtom object of the given atom container object
     *     as a property.
     *     <br>Then, all matches of the molecule's substructures with the internal sugar structure patterns are determined
     *     and in the following treated as candidate sugar structures. The candidates are then scanned for overlapping
     *     matches and the specific overlapping part (NOT the whole candidate) is removed from one of the overlapping
     *     pair. After this, there are no overlapping candidate sugar structures.
     *     <br>If the respective option is set, all candidates are also searched for substructures that are part of rings
     *     in the original molecule and these candidate structures (the WHOLE structure, not only the cyclic part) are
     *     discarded from the candidates.
     *     <br>If only terminal sugar moieties are to be removed, the sugar candidates are one-by-one tested for
     *     whether they are terminal or not and removed if they are. The procedure starts anew after iterating over all
     *     candidates and stops if no terminal sugar was removed in one whole iteration. The determination whether
     *     a sugar is terminal or not depends also on the set structure to keep mode because if the removal of a sugar
     *     candidate leads to an unconnected fragment the sugar candidate can still be terminal if the fragment is
     *     too small to keep and therefore discarded after the removal of the sugar. Resulting again in one connected
     *     remaining structure.
     *     <br>If also non-terminal sugars should be removed, simply all detected candidates are removed from the given
     *     molecular structure. If too small/light structures are also removed according to the current settings, the
     *     resulting structure will be cleared of all unconnected fragments that are not big/heavy enough. Still, a
     *     molecule consisting of multiple unconnected structures may result.
     *     <br>If the respective option is set, a property will be added to the given atom container specifying whether
     *     it contains (linear) sugar moieties or not.
     * </p>
     *
     * @param aMolecule the molecule to remove linear sugar moieties from
     * @param aShouldBeCloned true, if the sugar moieties should not be removed from the given atom container but a clone
     *                        of it should be generated and the sugars be removed from that
     * @return if the given atom container should NOT be cloned, this method returns the same given atom container after the
     * sugar removal; the returned molecule may be unconnected if also non-terminal sugars are removed according to
     * the settings and it may be empty if the resulting structure after sugar removal was too small to keep due to the
     * set structure to keep mode and the associated threshold (i.e. the molecule basically was a sugar)
     * @throws NullPointerException if the given atom container is 'null'
     * @throws CloneNotSupportedException if the given atom container should be cloned but does not allow cloning
     * @throws IllegalArgumentException if only terminal sugars should be removed but the given atom container already
     * contains multiple, unconnected structures which makes the determination of terminal and non-terminal structures
     * impossible
     */
    public IAtomContainer removeLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures. This makes the determination" +
                        "of terminal and non-terminal sugar moieties and therefore the sugar removal under the given " +
                        "settings impossible!");
            }
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        tmpNewMolecule = this.setIndices(tmpNewMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getLinearSugarCandidates(tmpNewMolecule);
        /*note: this means that there are matches of the linear sugar patterns and that they adhere to most of
        the given settings. The exception is that they might not be terminal*/
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsSugar);
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsSugar);
        }
        if (tmpContainsSugar) {
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpSugarCandidates);
            tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        }
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    //TODO: Try again to remove linear and circular candidates together (for a more valid terminal-check)
    // if it does not work, maybe chain multiple removals of linear and circular
    /**
     * Removes linear AND circular sugar moieties from the given atom container. Which substructures are removed depends
     * on the various settings available in this class.
     * <br>For the detailed procedure how the moieties are removed, see the respective methods for only linear and only
     * circular sugar removal. The detection of the two groups is separated but not their removal.
     * <br>If the respective option is set, a property will be added to the given atom container specifying whether
     * it contains (circular/linear) sugar moieties or not.
     *
     * @param aMolecule the molecule to remove sugar moieties from
     * @param aShouldBeCloned true, if the sugar moieties should not be removed from the given atom container but a clone
     *                        of it should be generated and the sugars be removed from that
     * @return if the given atom container should NOT be cloned, this method returns the same given atom container after the
     * sugar removal; the returned molecule may be unconnected if also non-terminal sugars are removed according to
     * the settings and it may be empty if the resulting structure after sugar removal was too small to keep due to the
     * set structure to keep mode and the associated threshold (i.e. the molecule basically was a sugar)
     * @throws NullPointerException if the given atom container is 'null'
     * @throws CloneNotSupportedException if the given atom container should be cloned but does not allow cloning
     * @throws IllegalArgumentException if only terminal sugars should be removed but the given atom container already
     * contains multiple, unconnected structures which makes the determination of terminal and non-terminal structures
     * impossible
     */
    public IAtomContainer removeAllSugars(IAtomContainer aMolecule, boolean aShouldBeCloned)
            throws NullPointerException, CloneNotSupportedException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        if (this.removeOnlyTerminal) {
            boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
            if (!tmpIsConnected) {
                throw new IllegalArgumentException("Only terminal sugar moieties should be removed but the given atom" +
                        "container already contains multiple unconnected structures. This makes the determination" +
                        "of terminal and non-terminal sugar moieties and therefore the sugar removal under the given " +
                        "settings impossible!");
            }
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        tmpNewMolecule = this.setIndices(tmpNewMolecule);
        //note: this has to be done stepwise because linear and circular sugar candidates can overlap
        List<IAtomContainer> tmpCircularSugarCandidates = this.getCircularSugarCandidates(tmpNewMolecule);
        boolean tmpContainsCircularSugars = !tmpCircularSugarCandidates.isEmpty();
        if (tmpContainsCircularSugars) {
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpCircularSugarCandidates);
            tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        }
        //exit here if molecule is empty after removal
        if (tmpNewMolecule.isEmpty()) {
            if (this.setPropertyOfSugarContainingMolecules) {
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsCircularSugars);
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsCircularSugars);
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, false);
            }
            return tmpNewMolecule;
        }
        //note: if only terminal sugars are removed, the atom container should not be disconnected at this point
        List<IAtomContainer> tmpLinearSugarCandidates = this.getLinearSugarCandidates(tmpNewMolecule);
        boolean tmpContainsLinearSugars = !tmpLinearSugarCandidates.isEmpty();
        if (tmpContainsLinearSugars) {
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpLinearSugarCandidates);
            tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        }
        boolean tmpContainsAnyTypeOfSugars = (tmpContainsCircularSugars || tmpContainsLinearSugars);
        if (this.setPropertyOfSugarContainingMolecules) {
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, tmpContainsAnyTypeOfSugars);
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, tmpContainsCircularSugars);
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, tmpContainsLinearSugars);
        }
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Utility method that can used to select the 'biggest' (i.e. the one with the highest heavy atom count) structure
     * from an atom container containing multiple unconnected structures, e.g. after the removal of both terminal and
     * non-terminal sugar moieties.
     * <br>The properties of the given atom container (IAtomContainer.getProperties()) are transferred to the returned
     * atom container.
     * <br>Note: This method does not clear away structures that are too small/light. It is independent of the settings
     * done in an object of this class.
     *
     * @param aMolecule the molecule to select the biggest structure from out of multiple unconnected structures
     * @return the biggest structure
     * @throws NullPointerException if the given atom container is 'null' or the CDK ConnectivityChecker is unable
     * to determine the unconnected structures
     */
    public static IAtomContainer selectBiggestUnconnectedFragment(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        if (tmpIsConnected) {
            return aMolecule;
        }
        Map<Object, Object> tmpProperties = aMolecule.getProperties();
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpBiggestFragment;
        if(tmpUnconnectedFragments != null && tmpUnconnectedFragments.getAtomContainerCount() >= 1) {
            tmpBiggestFragment = tmpUnconnectedFragments.getAtomContainer(0);
            int tmpBiggestFragmentHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(tmpBiggestFragment).size();
            for(IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()){
                int tmpFragmentHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(tmpFragment).size();
                if(tmpFragmentHeavyAtomCount > tmpBiggestFragmentHeavyAtomCount){
                    tmpBiggestFragment = tmpFragment;
                    tmpBiggestFragmentHeavyAtomCount = tmpFragmentHeavyAtomCount;
                }
            }
        } else {
            throw new NullPointerException("Could not detect the unconnected structures of the given atom container.");
        }
        tmpBiggestFragment.setProperties(tmpProperties);
        return tmpBiggestFragment;
    }

    /**
     * Utility method that can used to select the 'heaviest' (i.e. the one with the highest molecular weight) structure
     * from an atom container containing multiple unconnected structures, e.g. after the removal of both terminal and
     * non-terminal sugar moieties.
     * <br>The properties of the given atom container (IAtomContainer.getProperties()) are transferred to the returned
     * atom container.
     * <br>Note: This method does not clear away structures that are too small/light. It is independent of the settings
     * done in an object of this class.
     *
     * @param aMolecule the molecule to select the heaviest structure from out of multiple unconnected structures
     * @return the heaviest structure
     * @throws NullPointerException if the given atom container is 'null' or the CDK ConnectivityChecker is unable
     * to determine the unconnected structures
     */
    public static IAtomContainer selectHeaviestUnconnectedFragment(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        if (tmpIsConnected) {
            return aMolecule;
        }
        Map<Object, Object> tmpProperties = aMolecule.getProperties();
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpHeaviestFragment;
        if(tmpUnconnectedFragments != null && tmpUnconnectedFragments.getAtomContainerCount() >= 1) {
            tmpHeaviestFragment = tmpUnconnectedFragments.getAtomContainer(0);
            double tmpHeaviestFragmentWeight = AtomContainerManipulator.getMass(tmpHeaviestFragment);
            for(IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()){
                double tmpFragmentWeight = AtomContainerManipulator.getMass(tmpFragment);
                if(tmpFragmentWeight > tmpHeaviestFragmentWeight){
                    tmpHeaviestFragment = tmpFragment;
                    tmpHeaviestFragmentWeight = tmpFragmentWeight;
                }
            }
        } else {
            //if something went wrong
            return null;
        }
        tmpHeaviestFragment.setProperties(tmpProperties);
        return tmpHeaviestFragment;
    }

    /**
     * Utility method that can used to partition the unconnected structures in an atom container, e.g. after the removal
     * of both terminal and non-terminal sugar moieties, into a list of separate atom container objects and sort this
     * list with the following criteria with decreasing priority: atom count, molecular weight, bond count and sum of
     * bond orders.
     *
     * @param aMolecule the molecule whose unconnected structures to separate and sort
     * @return list of sorted atom containers representing the unconnected structures of the given molecule
     * @throws NullPointerException if the given atom container is 'null'
     */
    public static List<IAtomContainer> partitionAndSortUnconnectedFragments(IAtomContainer aMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        boolean tmpIsEmpty = aMolecule.isEmpty();
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        if (tmpIsConnected || tmpIsEmpty) {
            ArrayList<IAtomContainer> tmpFragmentList = new ArrayList<>(1);
            tmpFragmentList.add(aMolecule);
            return tmpFragmentList;
        }
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        int tmpSize = tmpUnconnectedFragments.getAtomContainerCount();
        ArrayList<IAtomContainer> tmpSortedList = new ArrayList<>(tmpSize);
        for (IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()) {
            tmpSortedList.add(tmpFragment);
        }
        /*Compares two IAtomContainers for order with the following criteria with decreasing priority:
            Compare atom count
            Compare molecular weight (heavy atoms only)
            Compare bond count
            Compare sum of bond orders (heavy atoms only)
        If no difference can be found with the above criteria, the IAtomContainers are considered equal.*/
        AtomContainerComparator tmpComparator = new AtomContainerComparator();
        tmpSortedList.sort(tmpComparator);
        return tmpSortedList;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Adds indices as properties to all atom objects of the given atom container to identify them uniquely within the
     * atom container and its clones. This is required e.g. for the determination of terminal vs. non-terminal sugar
     * moieties.
     *
     * @param aMolecule the molecule to process
     * @return the same atom container object
     * @throws NullPointerException if molecule is 'null' (note: no further parameter tests are implemented!)
     */
    private IAtomContainer setIndices(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int tmpIndex = 0;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpAtom.setProperty(SugarRemovalUtility.INDEX_PROPERTY_KEY, tmpIndex);
            tmpIndex++;
        }
        return aMolecule;
    }

    /**
     * Checks whether all exocyclic bonds connected to a given ring fragment of an original atom container are of single
     * order. The method iterates over all cyclic atoms and all of their bonds. So the runtime scales linear with the number
     * of cyclic atoms and their connected bonds.
     *
     * @param aRingToTest the ring fragment to test; exocyclic bonds do not have to be included in the fragment but if it
     *                    is a fused system of multiple rings, the internal interconnecting bonds of the different rings
     *                    need to be included; all its atoms need to be exactly the same objects as in the second atom
     *                    container parameter
     * @param anOriginalMolecule the molecule that contains the ring under investigation; The exocyclic bonds will be
     *                           queried from it
     * @return true, if all exocyclic bonds connected to the ring are of single order
     * @throws NullPointerException if a parameter is 'null' (note: no further parameter tests are implemented!)
     */
    private boolean areAllExocyclicBondsSingle(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule " +
                "is 'null'");
        int tmpAtomCountInRing = aRingToTest.getAtomCount();
        int tmpArrayListInitCapacity = tmpAtomCountInRing * 2;
        List<IBond> tmpExocyclicBondsList = new ArrayList<>(tmpArrayListInitCapacity);
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            List<IBond> tmpConnectedBondsList = anOriginalMolecule.getConnectedBondsList(tmpRingAtom);
            for (IBond tmpBond : tmpConnectedBondsList) {
                boolean tmpIsInRing = aRingToTest.contains(tmpBond);
                if (!tmpIsInRing) {
                    tmpExocyclicBondsList.add(tmpBond);
                }
            }
        }
        return (BondManipulator.getMaximumBondOrder(tmpExocyclicBondsList) == IBond.Order.SINGLE);
    }

    /**
     * Checks all exocyclic connections of the given ring to detect a glycosidic bond. Checklist for glycosidic bond:
     * Connected oxygen atom that is not in the ring, has two bonds that are both of single order and no bond partner
     * is a hydrogen atom.
     * <br>Note: This algorithm also classifies ester and ether bonds as a glycosidic bond which is not really bad. They
     * occur frequently in natural products.
     * <br>Note 2: The 'ring' is not tested for whether it is circular or not. So theoretically, this method
     * can also be used to detect glycosidic bonds of linear structures. BUT: The oxygen atom must not be part of the
     * structure itself. Due to the processing of candidate linear sugar moieties this can make it difficult to use
     * this method also for linear sugars.
     *
     * @param aRingToTest the candidate sugar ring
     * @param anOriginalMolecule the molecule in which the ring is contained as a substructure to query the connected
     *                           atoms from
     * @return true, if a glycosidic bond is detected
     * @throws NullPointerException if any parameter is 'null' (note: no further parameter tests are implemented!)
     */
    private boolean hasGlycosidicBond(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule)
            throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule " +
                "is 'null'");
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        boolean tmpContainsGlycosidicBond = false;
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            boolean tmpBreakOuterLoop = false;
            //check to avoid exceptions
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            List<IAtom> connectedAtomsList = anOriginalMolecule.getConnectedAtomsList(tmpRingAtom);
            for (IAtom tmpAtom : connectedAtomsList) {
                boolean tmpIsInRing = aRingToTest.contains(tmpAtom);
                if (!tmpIsInRing) {
                    String tmpSymbol = tmpAtom.getSymbol();
                    boolean tmpIsOxygen = (tmpSymbol.equals("O"));
                    if (tmpIsOxygen) {
                        List<IBond> tmpConnectedBondsList = anOriginalMolecule.getConnectedBondsList(tmpAtom);
                        boolean tmpHasOnlyTwoBonds = (tmpConnectedBondsList.size() == 2);
                        boolean tmpAllBondsAreSingle =
                                (BondManipulator.getMaximumBondOrder(tmpConnectedBondsList) == IBond.Order.SINGLE);
                        boolean tmpOneBondAtomIsHydrogen = false;
                        for (IBond tmpBond : tmpConnectedBondsList) {
                            for (IAtom tmpBondAtom : tmpBond.atoms()) {
                                if (tmpBondAtom.getSymbol().equals("H")) {
                                    tmpOneBondAtomIsHydrogen = true;
                                }
                            }
                        }
                        if ((tmpHasOnlyTwoBonds && tmpAllBondsAreSingle) && !tmpOneBondAtomIsHydrogen) {
                            tmpContainsGlycosidicBond = true;
                            tmpBreakOuterLoop = true;
                            break;
                        }
                    }
                }
            }
            if (tmpBreakOuterLoop) {
                break;
            }
        }
        return tmpContainsGlycosidicBond;
    }

    /**
     * Gives the number of attached exocyclic oxygen atoms of a given ring fragment of an original atom container.
     * The method iterates over all cyclic atoms and all of their connected atoms. So the runtime scales linear with
     * the number of cyclic atoms and their connected atoms.
     * <br>Note: The oxygen atoms are not tested for being attached by a single bond since in the algorithm, the
     * determination whether a candidate sugar ring has only exocyclic single bonds precedes the calling of this method.
     * <br>Note: The circularity of the given 'ring' is not tested, so this method could in theory also be used for linear
     * structures. But his does not make much sense.
     * <br>Note: This method does NOT check for hydroxy groups but for oxygen atoms. So e.g. the oxygen atom in a
     * glycosidic bond is counted.
     *
     * @param aRingToTest the ring fragment to test; exocyclic bonds do not have to be included in the fragment but if it
     *                    is a fused system of multiple rings, the internal interconnecting bonds of the different rings
     *                    need to be included; all its atoms need to be exactly the same objects as in the second atom
     *                    container parameter
     * @param anOriginalMolecule the molecule that contains the ring under investigation; The exocyclic bonds will be
     *                           queried from it
     * @return number of attached exocyclic oxygen atoms
     * @throws NullPointerException if a parameter is 'null' (note: no further parameter tests are implemented!)
     */
   private int getAttachedOxygenAtomCount(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule)
           throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule " +
                "is 'null'");
        int tmpExocyclicOxygenCounter = 0;
        Iterable<IAtom> tmpRingAtoms = aRingToTest.atoms();
        for (IAtom tmpRingAtom : tmpRingAtoms) {
            //check to avoid exceptions
            if (!anOriginalMolecule.contains(tmpRingAtom)) {
                continue;
            }
            List<IAtom> tmpConnectedAtomsList = anOriginalMolecule.getConnectedAtomsList(tmpRingAtom);
            for (IAtom tmpConnectedAtom : tmpConnectedAtomsList) {
                String tmpSymbol = tmpConnectedAtom.getSymbol();
                boolean tmpIsOxygen = tmpSymbol.equals("O");
                boolean tmpIsInRing = aRingToTest.contains(tmpConnectedAtom);
                if (tmpIsOxygen && !tmpIsInRing) {
                    tmpExocyclicOxygenCounter++;
                }
            }
        }
        return tmpExocyclicOxygenCounter;
    }

    /**
     * Simple decision making function for deciding whether a candidate sugar ring has enough attached, single-bonded
     * exocyclic oxygen atoms, called if this option is enabled in the current settings. The given number of oxygen atoms
     * is divided by the given number of atoms in the ring (should also contain the usually present oxygen atom in
     * a sugar ring) and the resulting ratio is checked for being equal or higher than the currently set
     * threshold.
     *
     * @param aNumberOfAtomsInRing number of atoms in the possible sugar ring, including the cyclic oxygen atom
     * @param aNumberOfAttachedExocyclicOxygenAtoms number of attached exocyclic oxygen atoms of the ring under
     *                                              investigation (if zero, false is returned)
     * @return true, if the calculated ratio is equal to or higher than the currently set threshold
     */
    private boolean doesRingHaveEnoughOxygenAtomsAttached(int aNumberOfAtomsInRing,
                                                          int aNumberOfAttachedExocyclicOxygenAtoms) {
        if (aNumberOfAtomsInRing == 0) {
            //better than throwing an exception here?
            return false;
        }
        double tmpAttachedOxygensToAtomsInRingRatio =
                ((double) aNumberOfAttachedExocyclicOxygenAtoms / (double) aNumberOfAtomsInRing);
        boolean tmpMeetsThreshold =
                (tmpAttachedOxygensToAtomsInRingRatio >= this.attachedOxygensToAtomsInRingRatioThreshold);
        return tmpMeetsThreshold;
    }

    /**
     * Removes all unconnected fragments that are too small/light according to the current structure to keep mode and
     * threshold setting. If all structures are too small/light, an empty atom container is returned. This does not
     * guarantee that the resulting atom container consists of only one connected structure.
     *
     * @param aMolecule the molecule to clear
     * @return the same atom container object, may be empty
     * @throws NullPointerException if molecule is 'null' (note: no further parameter tests are implemented!)
     */
    private IAtomContainer clearTooSmallStructures(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            return aMolecule;
        }
        IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        for (int i = 0; i < tmpComponents.getAtomContainerCount(); i++) {
            IAtomContainer tmpComponent = tmpComponents.getAtomContainer(i);
            boolean tmpIsTooSmall = this.isTooSmall(tmpComponent);
            if (tmpIsTooSmall) {
                //note: careful with removing things from sets/lists while iterating over it! But here it is ok because elements
                // are not removed from the same set that is iterated
                for (IAtom tmpAtom : tmpComponent.atoms()) {
                    //check to avoid exceptions
                    if (aMolecule.contains(tmpAtom)) {
                        aMolecule.removeAtom(tmpAtom);
                    }
                }
            }
        }
        //does not return null, but the atom container might be empty!
        return aMolecule;
    }

    /**
     * Checks whether the given molecule or structure is too small/light according to the current structure to keep mode
     * and threshold setting.
     *
     * @param aMolecule the molecule to check
     * @return true, if the given structure is too light/small
     * @throws NullPointerException if molecule is 'null' (note: no further parameter tests are implemented!)
     * @throws UnsupportedOperationException if an unknown StructureToKeepMode enum constant is set
     */
    private boolean isTooSmall(IAtomContainer aMolecule) throws NullPointerException, UnsupportedOperationException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        boolean tmpIsTooSmall;
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            tmpIsTooSmall = false;
        } else if (this.structuresToKeepMode == StructuresToKeepMode.HEAVY_ATOM_COUNT) {
            int tmpHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(aMolecule).size();
            tmpIsTooSmall = tmpHeavyAtomCount < this.structureToKeepModeThreshold;
        } else if (this.structuresToKeepMode == StructuresToKeepMode.MOLECULAR_WEIGHT) {
            double tmpMolWeight = AtomContainerManipulator.getMass(aMolecule, AtomContainerManipulator.MolWeight);
            tmpIsTooSmall = tmpMolWeight < this.structureToKeepModeThreshold;
        } else {
            throw new UnsupportedOperationException("Undefined StructuresToKeepMode setting!");
        }
        return tmpIsTooSmall;
    }

    /**
     * Checks whether the given substructure is terminal in the given parent molecule. To do this, the substructure and
     * the parent molecule are cloned, the substructure is removed in the parent molecule clone (identified via the index
     * properties) and finally it is checked whether the parent molecule clone still consists of only one connected
     * structure. If that is the case, the substructure is terminal. If the structure to keep mode is not set to 'keep
     * all structures', too small resulting fragments are cleared away in between.
     *
     * @param aSubstructure the substructure to check for whether it is terminal
     * @param aParentMolecule the molecule the substructure is a part of
     * @return true, if the substructure is terminal
     * @throws NullPointerException if any parameter is 'null'
     * @throws IllegalArgumentException if the substructure is not part of the parent molecule or if the parent molecule
     * is already unconnected (i.e. consists of multiple, unconnected substructures)
     * @throws CloneNotSupportedException if one of the atom containers cannot be cloned
     */
    private boolean isTerminal(IAtomContainer aSubstructure, IAtomContainer aParentMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aSubstructure, "Given substructure is 'null'.");
        Objects.requireNonNull(aParentMolecule, "Given parent molecule is 'null'.");
        boolean tmpIsParent = true;
        for (IAtom tmpAtom : aSubstructure.atoms()) {
            if (!aParentMolecule.contains(tmpAtom)) {
                tmpIsParent = false;
                break;
            }
        }
        if (!tmpIsParent) {
            throw new IllegalArgumentException("Given substructure is not part of the given parent molecule.");
        }
        boolean tmpIsUnconnected = !ConnectivityChecker.isConnected(aParentMolecule);
        if (tmpIsUnconnected) {
            throw new IllegalArgumentException("Parent molecule is already unconnected.");
        }
        //</editor-fold>
        //TODO/discuss: Is there a better way?
        boolean tmpIsTerminal;
        IAtomContainer tmpMoleculeClone = aParentMolecule.clone();
        IAtomContainer tmpSubstructureClone = aSubstructure.clone();
        HashMap<Integer, IAtom> tmpIndexToAtomMap = new HashMap<>(tmpMoleculeClone.getAtomCount() + 1, 1);
        for (IAtom tmpAtom : tmpMoleculeClone.atoms()) {
            tmpIndexToAtomMap.put(tmpAtom.getProperty(SugarRemovalUtility.INDEX_PROPERTY_KEY), tmpAtom);
        }
        for (IAtom tmpAtom : tmpSubstructureClone.atoms()) {
            tmpMoleculeClone.removeAtom(tmpIndexToAtomMap.get(tmpAtom.getProperty(SugarRemovalUtility.INDEX_PROPERTY_KEY)));
        }
        boolean tmpIsConnected = ConnectivityChecker.isConnected(tmpMoleculeClone);
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            tmpIsTerminal = tmpIsConnected;
        } else {
            if (tmpIsConnected) {
                tmpIsTerminal = true;
            } else {
                IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(tmpMoleculeClone);
                for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
                    boolean tmpIsTooSmall = this.isTooSmall(tmpComponent);
                    if (tmpIsTooSmall) {
                        tmpMoleculeClone.remove(tmpComponent);
                    }
                }
                tmpIsTerminal = ConnectivityChecker.isConnected(tmpMoleculeClone);
            }
        }
        return tmpIsTerminal;
    }

    /**
     * Removes the given sugar candidate structures (circular and linear) from the given molecule. The candidates have
     * to be pre-filtered according to the current settings. The only setting influencing the removal here is the
     * option to remove only terminal sugar moieties (which is internally influenced by the structure to keep mode and
     * threshold).
     * <br>If only terminal sugar moieties are to be removed, the sugar candidates are one-by-one tested for
     * whether they are terminal or not and removed if they are. The procedure starts anew after iterating over all
     * candidates and stops if no terminal sugar was removed in one whole iteration. Too small/light, unconnected
     * structures are cleared after every removal (but not at the end!).
     * <br>If also non-terminal sugars should be removed, simply all detected candidates are removed from the given
     * molecular structure.
     * <br>To clear away unconnected, too small/light structures and generate valid valences at the end, the
     * post-processing method must be called on the resulting atom container.
     *
     * @param aMolecule the molecule to remove the sugar candidates from
     * @param aCandidateList the list of candidate sugar moieties in the given molecule
     * @return the given molecule (i.e. atom container object) after removal of the sugar moieties
     * @throws NullPointerException if any parameter is 'null'
     * @throws IllegalArgumentException if at least one atom in the candidate list is not actually part of the molecule
     */
    private IAtomContainer removeSugarCandidates(IAtomContainer aMolecule, List<IAtomContainer> aCandidateList)
            throws NullPointerException, IllegalArgumentException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        Objects.requireNonNull(aCandidateList, "Given list is 'null'.");
        if (aCandidateList.isEmpty()) {
            return aMolecule;
        }
        for (IAtomContainer tmpSubstructure : aCandidateList) {
            boolean tmpIsParent = true;
            for (IAtom tmpAtom : tmpSubstructure.atoms()) {
                if (!aMolecule.contains(tmpAtom)) {
                    tmpIsParent = false;
                    break;
                }
            }
            if (!tmpIsParent) {
                throw new IllegalArgumentException("At least one of the possible sugar-like substructures is not " +
                        "actually part of the given molecule.");
            }
        }
        //</editor-fold>
        IAtomContainer tmpNewMolecule = aMolecule;
        List<IAtomContainer> tmpSugarCandidates = aCandidateList;
        if (this.removeOnlyTerminal) {
            //Only terminal sugars should be removed
            //but the definition of terminal depends on the set structures to keep mode!
            //decisions based on this setting are made in the respective private method
            //No unconnected structures result at the end or at an intermediate step
            boolean tmpContainsNoTerminalSugar = false;
            while (!tmpContainsNoTerminalSugar) {
                boolean tmpSomethingWasRemoved = false;
                for (int i = 0; i < tmpSugarCandidates.size(); i++) {
                    IAtomContainer tmpCandidate = tmpSugarCandidates.get(i);
                    boolean tmpIsTerminal = false;
                    try {
                        tmpIsTerminal = this.isTerminal(tmpCandidate, tmpNewMolecule);
                    } catch (CloneNotSupportedException aCloneNotSupportedException) {
                        SugarRemovalUtility.LOGGER.log(Level.WARNING, aCloneNotSupportedException.toString(),
                                aCloneNotSupportedException);
                        throw new IllegalArgumentException("Could not clone one candidate and therefore not determine " +
                                "whether it is terminal or not.");
                    }
                    if (tmpIsTerminal) {
                        for (IAtom tmpAtom : tmpCandidate.atoms()) {
                            if (tmpNewMolecule.contains(tmpAtom)) {
                                tmpNewMolecule.removeAtom(tmpAtom);
                            }
                        }
                        //TODO: Post-process after every removal? Formerly non-terminal sugars might have only attached
                        // oxygen atoms now instead of hydroxy groups. This might be a problem!
                        tmpSugarCandidates.remove(i);
                        //The removal shifts the remaining indices!
                        i = i - 1;
                        //to clear away leftover unconnected fragments that are not to be kept due to the settings
                        tmpNewMolecule = this.clearTooSmallStructures(tmpNewMolecule);
                        //atom container may be empty after that
                        if (tmpNewMolecule.isEmpty()) {
                            tmpContainsNoTerminalSugar = true;
                            break;
                        }
                        tmpSomethingWasRemoved = true;
                    }
                }
                if (!tmpSomethingWasRemoved) {
                    tmpContainsNoTerminalSugar = true;
                }
            }
        } else {
            //all sugar moieties are removed, may result in an unconnected atom container
            for (IAtomContainer tmpSugarCandidate : tmpSugarCandidates) {
                for (IAtom tmpAtom : tmpSugarCandidate.atoms()) {
                    if (tmpNewMolecule.contains(tmpAtom)) {
                        tmpNewMolecule.removeAtom(tmpAtom);
                    }
                }
            }
        }
        //To clear away unconnected, too small structures and generate valid valences, the post-processing method must
        // be called
        return tmpNewMolecule;
    }

    /**
     * Clears away too small/light structures from the given molecule (may result in an empty atom container) and
     * generates valid valences on the remaining molecule by the addition of implicit hydrogen atoms to open valences.
     *
     * @param aMolecule the molecule to post-process
     * @return the same molecule (i.e. atom container object)
     * @throws NullPointerException if the given molecule is 'null' (note: no further parameter tests are implemented!)
     */
    private IAtomContainer postProcessAfterRemoval(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        IAtomContainer tmpNewMolecule = aMolecule;
        //if too small / too light, unconnected structures should be discarded, this is done now
        //otherwise, the possibly unconnected atom container is returned
        //Even if only terminal sugars are removed, the resulting, connected structure may still be too small to keep!
        if (this.structuresToKeepMode != StructuresToKeepMode.ALL) {
            tmpNewMolecule = this.clearTooSmallStructures(tmpNewMolecule);
        }
        if (!tmpNewMolecule.isEmpty()) {
            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpNewMolecule);
                CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(tmpNewMolecule);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
            }
        }
        return tmpNewMolecule;
    }

    /**
     * Extracts circular sugar candidate structures from the given molecule according to the current settings.
     * <br>First, the isolated rings of the molecule are matched with the internal circular sugar structure patterns.
     * All matching rings are candidate sugar structures.
     * <br>Discarded from the candidates are rings that have an exocyclic double or triple bond. If the respective
     * option is set, all candidates are discarded that do not have a glycosidic bond. Furthermore, all
     * candidates are discarded that have too few attached exocyclic oxygen atoms if the respective option is set.
     * <br>The candidates returned here might not be terminal!
     *
     * @param aMolecule the molecule to extract circular sugar candidates from
     * @return a list of substructures in the given molecule that are regarded as circular sugar moieties
     * @throws NullPointerException if the given molecule is 'null' (note: no further parameter tests are implemented!)
     */
    private List<IAtomContainer> getCircularSugarCandidates(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        IAtomContainer tmpNewMolecule = aMolecule;
        int[][] tmpAdjList = GraphUtil.toAdjList(tmpNewMolecule);
        //efficient computation/partitioning of the ring systems
        RingSearch tmpRingSearch = new RingSearch(tmpNewMolecule, tmpAdjList);
        List<IAtomContainer> tmpIsolatedRings = tmpRingSearch.isolatedRingFragments();
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(tmpIsolatedRings.size());
        for(IAtomContainer tmpReferenceRing : this.ringSugars){
            for(IAtomContainer tmpIsolatedRing : tmpIsolatedRings){
                boolean tmpIsIsomorph = false;
                try {
                    tmpIsIsomorph = this.univIsomorphismTester.isIsomorph(tmpReferenceRing, tmpIsolatedRing);
                } catch (CDKException aCDKException) {
                    SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
                    continue;
                }
                if (tmpIsIsomorph) {
                    /*note: another requirement of a suspected sugar ring is that it contains only single bonds.
                     * This is not tested here because all the structures in the reference rings do meet this criterion.
                     * But a structure that does not meet this criterion could be added to the references by the user.*/
                    //do not remove rings with non-single exocyclic bonds, they are not sugars (not an option!)
                    boolean tmpAreAllExocyclicBondsSingle = this.areAllExocyclicBondsSingle(tmpIsolatedRing, tmpNewMolecule);
                    if (!tmpAreAllExocyclicBondsSingle) {
                        continue;
                    }
                    //do not remove rings without an attached glycosidic bond if this option is set
                    if (this.detectGlycosidicBond) {
                        boolean tmpHasGlycosidicBond = this.hasGlycosidicBond(tmpIsolatedRing, tmpNewMolecule);
                        if (!tmpHasGlycosidicBond) {
                            continue;
                        }
                    }
                    //do not remove rings with 'too few' attached oxygens if this option is set
                    if (this.includeNrOfAttachedOxygens) {
                        int tmpExocyclicOxygenCount = this.getAttachedOxygenAtomCount(tmpIsolatedRing, tmpNewMolecule);
                        int tmpAtomsInRing = tmpIsolatedRing.getAtomCount();
                        boolean tmpAreEnoughOxygensAttached = this.doesRingHaveEnoughOxygenAtomsAttached(tmpAtomsInRing,
                                tmpExocyclicOxygenCount);
                        if (!tmpAreEnoughOxygensAttached) {
                            continue;
                        }
                    }
                    //if sugar ring has not been excluded yet, the molecule contains sugars, although they might not
                    // be terminal
                    tmpSugarCandidates.add(tmpIsolatedRing);
                }
            }
        }
        return tmpSugarCandidates;
    }

    /**
     * Extracts linear sugar candidate structures from the given molecule according to the current settings.
     * <br>First, all matches of the molecule's substructures with the internal sugar structure patterns are determined
     * and in the following treated as candidate sugar structures. The candidates are then scanned for overlapping
     * matches and the specific overlapping part (NOT the whole candidate) is removed from one of the overlapping
     * pair. After this, there are no overlapping candidate sugar structures.
     * <br>If the respective option is set, all candidates are also searched for substructures that are part of rings
     * in the original molecule and these atoms (and bonds, NOT the whole candidate) are discarded from the candidates.
     * <br>The candidates returned here might not be terminal!
     *
     * @param aMolecule the molecule to extract linear sugar candidates from
     * @return a list of substructures in the given molecule that are regarded as linear sugar moieties
     * @throws NullPointerException if the given molecule is 'null' (note: no further parameter tests are implemented!)
     */
    private List<IAtomContainer> getLinearSugarCandidates(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        IAtomContainer tmpNewMolecule = aMolecule;
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(aMolecule.getAtomCount() / 2);
        for(DfPattern tmpLinearSugarPattern : this.linearSugarPatterns) {
            /*unique in this case means that the same match cannot be in this collection multiple times but they can
            still overlap! Overlapping atoms are removed in the following lines.*/
            Mappings tmpMappings = tmpLinearSugarPattern.matchAll(tmpNewMolecule);
            Mappings tmpUniqueMappings = tmpMappings.uniqueAtoms();
            Iterable<IAtomContainer> tmpUniqueSubstructureMappings = tmpUniqueMappings.toSubstructures();
            for (IAtomContainer tmpMatchedStructure : tmpUniqueSubstructureMappings) {
                tmpSugarCandidates.add(tmpMatchedStructure);
            }
        }
        //TODO/discuss: Is there a better way to get non-overlapping matches?
        HashSet<Integer> tmpSugarCandidateAtomsSet = new HashSet<>(aMolecule.getAtomCount() + 2, 1);
        for (int i = 0; i < tmpSugarCandidates.size(); i++) {
            IAtomContainer tmpCandidate = tmpSugarCandidates.get(i);
            for (int j = 0; j < tmpCandidate.getAtomCount(); j++) {
                IAtom tmpAtom = tmpCandidate.getAtom(j);
                int tmpAtomIndex = tmpAtom.getProperty(SugarRemovalUtility.INDEX_PROPERTY_KEY);
                boolean tmpIsAtomAlreadyInCandidates = tmpSugarCandidateAtomsSet.contains(tmpAtomIndex);
                if (tmpIsAtomAlreadyInCandidates) {
                    tmpCandidate.removeAtom(tmpAtom);
                    //The removal shifts the remaining indices!
                    j = j - 1;
                } else {
                    tmpSugarCandidateAtomsSet.add(tmpAtomIndex);
                }
            }
            if (tmpCandidate.isEmpty()) {
                tmpSugarCandidates.remove(tmpCandidate);
                //The removal shifts the remaining indices!
                i = i - 1;
            }
        }
        if (!this.removeLinearSugarsInRing && !tmpSugarCandidates.isEmpty()) {
            int[][] tmpAdjList = GraphUtil.toAdjList(tmpNewMolecule);
            RingSearch tmpRingSearch = new RingSearch(tmpNewMolecule, tmpAdjList);
            for (int i = 0; i < tmpSugarCandidates.size(); i++) {
                IAtomContainer tmpCandidate = tmpSugarCandidates.get(i);
                for (int j = 0; j < tmpCandidate.getAtomCount(); j++) {
                    IAtom tmpAtom = tmpCandidate.getAtom(j);
                    if (tmpRingSearch.cyclic(tmpAtom)) {
                        //  *this is the routine to remove the whole candidate*
                        tmpSugarCandidates.remove(i);
                        //removal shifts the remaining indices
                        i = i - 1;
                        break;
                        //  *this would be the routine to remove only the cyclic atoms (see also add. part below)*
                        /*if (tmpCandidate.contains(tmpAtom)) {
                            tmpCandidate.removeAtom(tmpAtom);
                            //The removal shifts the remaining indices!
                            j = j - 1;
                        }*/
                    }
                }
                //  *this would be the routine to remove only the cyclic atoms (part 2)*
                /*if (tmpCandidate.isEmpty()) {
                    tmpSugarCandidates.remove(i);
                    //The removal shifts the remaining indices!
                    i = i - 1;
                }*/
            }
        }
        return tmpSugarCandidates;
    }
    //</editor-fold>
}