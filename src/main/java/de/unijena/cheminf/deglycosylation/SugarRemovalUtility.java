/**
 * TODO: Add license header (switch to MIT?)
 */
package de.unijena.cheminf.deglycosylation;

/**
 * TODO:
 * - write docs
 *
 * To discuss:
 * - If all remaining structures are kept, hydroxy groups of the sugars are also kept. If only terminal sugars should be
 * removed, no sugar is removed in this case! So remove sugars with attached hydroxy groups also?
 * - Note: The situation at the previous connection point is unclear. When circular sugars are removed, a hydroxy
 * group remains at the core structure (if there was a glycosidic bond). But for the linear sugars, no general statement
 * like this can be made.
 * - add hetero atom count as StructureToKeepMode option?
 * - add detection of glycosidic bond for linear sugars? The respective method would already be suitable for linear sugars
 * - remove deprecated method?
 * - license
 * - make a console application to deploy as jar? Inputs and outputs?
 * - add more linear sugars?
 *      - e.g. aldopentose, different deoxypentoses/deoxyhexoses, heptoses, octoses (?)
 * - Do we need a general linear sugar pattern?
 *      - e.g. a single-bonded, simple carbon chain of sufficient length where nearly all carbon atoms have one hydroxy
 *      or keto group
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
import java.util.ConcurrentModificationException;
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
     * Property key for index that is added to any IAtom object in a given IAtomContainer object for unique
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
            //"OCC(O)C(O)C(O)C(O)CO", //hexitol TODO/discuss: this is a doublet!
            //"O=CC(O)C(O)C(O)C(O)CO", //aldohexose TODO/discuss: this is a doublet!
            "CCCCC(O)C(=O)O", //2-hydroxyhexanoic acid TODO/discuss: Is this a sugar?
            "CC(=O)CC(=O)CCC(=O)O", //4,6-dioxoheptanoic acid TODO/discuss: Is this a sugar?
            "O=C(O)CC(O)CC(=O)O", //3-hydroxypentanedioic acid TODO/discuss: Is this a sugar?
            "O=C(O)C(=O)C(=O)C(O)C(O)CO", //hexo-2,3-diulosonic acid
            "O=C(O)CCC(O)C(=O)O", //2-hydroxypentanedioic acid TODO/discuss: Is this a sugar?
            //"O=CC(O)C(O)C(O)C(O)CO", //aldohexose TODO/discuss: this is a doublet!
            "O=C(CO)C(O)C(O)CO" //ketopentose
    };

    /**
     * TODO: Add doc
     */
    public static final String [] RING_SUGARS_SMILES = {
            "C1CCOC1", //tetrahydrofuran to match all 5-membered sugar rings, formed by pentoses
            "C1CCOCC1", //tetrahydropyran to match all 6-membered sugar rings, formed by hexoses
            "C1CCCOCC1" //oxepane to match all 7-membered sugar rings, formed by heptoses
    };

    /**
     * TODO
     */
    public static final boolean DETECT_GLYCOSIDIC_BOND_DEFAULT = false;

    /**
     * TODO
     */
    public static final boolean REMOVE_ONLY_TERMINAL_DEFAULT = true;

    /**
     * TODO
     */
    public static final StructuresToKeepMode STRUCTURES_TO_KEEP_MODE_DEFAULT = StructuresToKeepMode.HEAVY_ATOM_COUNT;

    /**
     * TODO
     */
    public static final boolean INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT = true;

    /**
     * TODO
     */
    public static final double ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT = 0.5;

    /**
     * TODO
     */
    public static final boolean REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT = false;

    /**
     * TODO
     */
    public static final boolean SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT = true;
    //</editor-fold>
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(SugarRemovalUtility.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * TODO: Add doc
     */
    private final UniversalIsomorphismTester univIsomorphismTester;
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    /**
     * TODO: Add doc
     */
    private List<IAtomContainer> linearSugars;

    /**
     * TODO: Add doc
     */
    private List<IAtomContainer> ringSugars;

    /**
     * TODO: Add doc
     */
    private List<DfPattern> linearSugarPatterns;

    /**
     * TODO: Add doc
     */
    private boolean detectGlycosidicBond;

    /**
     * TODO
     */
    private boolean removeOnlyTerminal;

    /**
     * TODO
     */
    private StructuresToKeepMode structuresToKeepMode;

    /**
     * TODO
     */
    private int mwOrHacThreshold;

    /**
     * TODO
     */
    private boolean includeNrOfAttachedOxygens;

    /**
     * TODO
     */
    private double attachedOxygensToAtomsInRingRatioThreshold;

    /**
     * TODO
     */
    private boolean removeLinearSugarsInRing;

    /**
     * TODO
     */
    private boolean setPropertyOfSugarContainingMolecules;
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * TODO: Add doc
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
        this.mwOrHacThreshold = this.structuresToKeepMode.defaultThreshold;
        this.includeNrOfAttachedOxygens = SugarRemovalUtility.INCLUDE_NR_OF_ATTACHED_OXYGEN_DEFAULT;
        this.attachedOxygensToAtomsInRingRatioThreshold = SugarRemovalUtility.ATTACHED_OXYGENS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT;
        this.removeLinearSugarsInRing = SugarRemovalUtility.REMOVE_LINEAR_SUGARS_IN_RING_DEFAULT;
        this.setPropertyOfSugarContainingMolecules = SugarRemovalUtility.SET_PROPERTY_OF_SUGAR_CONTAINING_MOLECULES_DEFAULT;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get/is">
    /**
     * TODO
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
            tmpSmilesList.add(tmpSmiles);
        }
        return tmpSmilesList;
    }

    /**
     * TODO
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
            tmpSmilesList.add(tmpSmiles);
        }
        return tmpSmilesList;
    }

    /**
     * TODO
     * @return
     */
    public boolean isGlycosidicBondDetected() {
        return this.detectGlycosidicBond;
    }

    /**
     * TODO
     */
    public boolean areOnlyTerminalSugarsRemoved() {
        return this.removeOnlyTerminal;
    }

    /**
     * TODO
     */
    public StructuresToKeepMode getStructuresToKeepMode() {
        return this.structuresToKeepMode;
    }

    /**
     * TODO
     */
    public int getMwOrHacThreshold() {
        return this.mwOrHacThreshold;
    }

    /**
     * TODO
     */
    public boolean isNrOfAttachedOxygensIncluded() {
        return this.includeNrOfAttachedOxygens;
    }

    /**
     * TODO
     */
    public double getAttachedOxygensToAtomsInRingRatioThreshold() {
        return this.attachedOxygensToAtomsInRingRatioThreshold;
    }

    /**
     * TODO
     */
    public boolean areLinearSugarsInRingsRemoved() {
        return this.removeLinearSugarsInRing;
    }

    /**
     * TODO
     */
    public boolean isPropertyOfSugarContainingMoleculesSet() {
        return this.setPropertyOfSugarContainingMolecules;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set/add/clear">
    /**
     * TODO: Add doc
     */
    public void addCircularSugar(IAtomContainer aCircularSugar) throws NullPointerException, IllegalArgumentException, IllegalStateException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aCircularSugar, "Given atom container is 'null'");
        if (aCircularSugar.isEmpty()) {
            throw new IllegalArgumentException("Given atom container is empty.");
        }
        boolean tmpSugarAlreadyPresent = false;
        UniversalIsomorphismTester tmpUnivIsomorphTester = new UniversalIsomorphismTester();
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
     * TODO
     */
    public void addCircularSugar(String aSmilesCode) throws NullPointerException, CDKException, IllegalArgumentException, IllegalStateException {
        Objects.requireNonNull(aSmilesCode, "Given SMILES code is 'null'");
        if (aSmilesCode.isEmpty()) {
            throw new IllegalArgumentException("Given SMILES code is empty");
        }
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpRingSugar = tmpSmiPar.parseSmiles(aSmilesCode);
        this.addCircularSugar(tmpRingSugar);
    }

    /**
     * TODO
     */
    public void addLinearSugar(IAtomContainer aLinearSugar) throws NullPointerException, IllegalArgumentException, IllegalStateException {
        //<editor-fold desc="Checks">
        Objects.requireNonNull(aLinearSugar, "Given atom container is 'null'");
        if (aLinearSugar.isEmpty()) {
            throw new IllegalArgumentException("Given atom container is empty.");
        }
        boolean tmpSugarAlreadyPresent = false;
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
     * TODO
     */
    public void addLinearSugar(String aSmilesCode) throws NullPointerException, CDKException, IllegalArgumentException, IllegalStateException {
        Objects.requireNonNull(aSmilesCode, "Given SMILES code is 'null'");
        if (aSmilesCode.isEmpty()) {
            throw new IllegalArgumentException("Given SMILES code is empty");
        }
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpLinearSugar = tmpSmiPar.parseSmiles(aSmilesCode);
        this.addLinearSugar(tmpLinearSugar);
    }

    /**
     * TODO
     */
    public void clearCircularSugars() {
        this.ringSugars.clear();
    }

    /**
     * TODO
     */
    public void clearLinearSugars() {
        this.linearSugars.clear();
        this.linearSugarPatterns.clear();
    }

    /**
     * TODO
     */
    public void setDetectGlycosidicBond(boolean aBoolean) {
        this.detectGlycosidicBond = aBoolean;
    }

    /**
     * TODO
     */
    public void setRemoveOnlyTerminalSugars(boolean aBoolean) {
        this.removeOnlyTerminal = aBoolean;
    }

    /**
     * TODO, add note that it overrides the previously set threshold!!
     */
    public void setStructuresToKeepMode(StructuresToKeepMode aMode) throws NullPointerException {
        Objects.requireNonNull(aMode, "Given mode is 'null'.");
        this.structuresToKeepMode = aMode;
        this.mwOrHacThreshold = this.structuresToKeepMode.getDefaultThreshold();
    }

    /**
     * TODO
     */
    public void setStructuresToKeepThreshold(int aThreshold) throws IllegalArgumentException {
        if ((this.structuresToKeepMode == StructuresToKeepMode.ALL) && (aThreshold > 1)) {
            throw new IllegalArgumentException("The mode is currently set to keep all structures, so a threshold > 1 makes no sense.");
        }
        if (aThreshold < 0) {
            throw new IllegalArgumentException("Threshold cannot be negative.");
        }
        this.mwOrHacThreshold = aThreshold;
    }

    /**
     * TODO
     */
    public void setIncludeNrOfAttachedOxygens(boolean aBoolean) {
        this.includeNrOfAttachedOxygens = aBoolean;
    }

    /**
     * TODO
     */
    public void setAttachedOxygensToAtomsInRingRatioThreshold(double aDouble) throws IllegalArgumentException {
        //<editor-fold desc="Checks">
        boolean tmpIsFinite = Double.isFinite(aDouble); //false for NaN and infinity arguments
        boolean tmpIsNegative = (aDouble < 0);
        if(!tmpIsFinite || tmpIsNegative) {
            throw new IllegalArgumentException("Given double is NaN, infinite or negative.");
        }
        boolean tmpIsBiggerThanOne = (aDouble > 1);
        if (tmpIsBiggerThanOne) {
            throw new IllegalArgumentException("Threshold cannot be bigger than one.");
        }
        if ((this.includeNrOfAttachedOxygens == false) && aDouble > 0) {
            throw new IllegalArgumentException("The number of attached oxygen atoms is currently not included in the " +
                    "decision making process, so a ration threshold > 0 makes no sense.");
        }
        //</editor-fold>
        this.attachedOxygensToAtomsInRingRatioThreshold = aDouble;
    }

    /**
     * TODO
     */
    public void setRemoveLinearSugarsInRing(boolean aBoolean) {
        this.removeLinearSugarsInRing = aBoolean;
    }

    /**
     * TODO
     */
    public void setPropertyOfSugarContainingMolecules(boolean aBoolean) {
        this.setPropertyOfSugarContainingMolecules = aBoolean;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * TODO
     */
    public boolean hasLinearSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        aMolecule = this.setIndices(aMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getLinearSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (tmpContainsSugar && this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, true);
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, true);
        }
        return tmpContainsSugar;
    }

    /**
     * TODO
     */
    public boolean hasCircularSugars(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return false;
        }
        aMolecule = this.setIndices(aMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (tmpContainsSugar && this.setPropertyOfSugarContainingMolecules) {
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, true);
            aMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, true);
        }
        return tmpContainsSugar;
    }

    /**
     * TODO
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
        if (tmpContainsSugar && this.setPropertyOfSugarContainingMolecules) {
            if (tmpContainsCircularSugar || tmpContainsLinearSugar) {
                aMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, true);
            }
            if (tmpContainsCircularSugar) {
                aMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, true);
            }
            if (tmpContainsLinearSugar) {
                aMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, true);
            }
        }
        return tmpContainsSugar;
    }

    /**
     * TODO
     */
    public IAtomContainer removeCircularSugars(IAtomContainer aMolecule, boolean aShouldBeCloned) throws NullPointerException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        tmpNewMolecule = this.setIndices(tmpNewMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (tmpContainsSugar) {
            /*note: this means that there are matches of the circular sugar patterns and that they adhere to most of
            the given settings. The exception is that they might not be terminal*/
            if (this.setPropertyOfSugarContainingMolecules) {
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, true);
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, true);
            }
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpSugarCandidates);
            tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        }
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    /**
     * TODO
     */
    public IAtomContainer removeLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned) throws NullPointerException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        tmpNewMolecule = this.setIndices(tmpNewMolecule);
        List<IAtomContainer> tmpSugarCandidates = this.getLinearSugarCandidates(tmpNewMolecule);
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (tmpContainsSugar) {
            /*note: this means that there are matches of the linear sugar patterns and that they adhere to most of
            the given settings. The exception is that they might not be terminal*/
            if (this.setPropertyOfSugarContainingMolecules) {
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, true);
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, true);
            }
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpSugarCandidates);
            tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        }
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    /**
     * TODO
     */
    public IAtomContainer removeAllSugars(IAtomContainer aMolecule, boolean aShouldBeCloned) throws NullPointerException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (aMolecule.isEmpty()) {
            return aMolecule;
        }
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        tmpNewMolecule = this.setIndices(tmpNewMolecule);
        List<IAtomContainer> tmpCircularSugarCandidates = this.getCircularSugarCandidates(aMolecule);
        boolean tmpContainsCircularSugars = !tmpCircularSugarCandidates.isEmpty();
        List<IAtomContainer> tmpLinearSugarCandidates = this.getLinearSugarCandidates(tmpNewMolecule);
        boolean tmpContainsLinearSugars = !tmpLinearSugarCandidates.isEmpty();
        if (this.setPropertyOfSugarContainingMolecules) {
            if (tmpContainsCircularSugars || tmpContainsLinearSugars) {
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY, true);
            }
            if (tmpContainsCircularSugars) {
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, true);
            }
            if (tmpContainsLinearSugars) {
                tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, true);
            }
        }
        if (tmpContainsCircularSugars && tmpContainsLinearSugars) {
            List<IAtomContainer> tmpAllSugarCandidates = new ArrayList<>(tmpCircularSugarCandidates.size()
                    + tmpLinearSugarCandidates.size());
            boolean tmpAddingCircularsWasSuccessful = tmpAllSugarCandidates.addAll(tmpCircularSugarCandidates);
            boolean tmpAddingLinearsWasSuccessful = tmpAllSugarCandidates.addAll(tmpLinearSugarCandidates);
            if (tmpAddingCircularsWasSuccessful && tmpAddingLinearsWasSuccessful) {
                tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpAllSugarCandidates);
                tmpContainsCircularSugars = false;
                tmpContainsLinearSugars = false;
            }
        }
        if (tmpContainsCircularSugars) {
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpCircularSugarCandidates);
        }
        if (tmpContainsLinearSugars) {
            tmpNewMolecule = this.removeSugarCandidates(tmpNewMolecule, tmpLinearSugarCandidates);
        }
        tmpNewMolecule = this.postProcessAfterRemoval(tmpNewMolecule);
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    /**
     * TODO: Add doc
     *
     * @param aMolecule
     * @return
     * @throws
     */
    @Deprecated
    public IAtomContainer removeSugars(IAtomContainer aMolecule, boolean aShouldBeCloned) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        IAtomContainer tmpNewMolecule;
        try {
            if (aShouldBeCloned) {
                tmpNewMolecule = aMolecule.clone();
            } else {
                tmpNewMolecule = aMolecule;
            }
            //***removing ring sugars***
            int[][] tmpAdjList = GraphUtil.toAdjList(tmpNewMolecule);
            //efficient computation/partitioning of the ring systems
            RingSearch tmpRingSearch = new RingSearch(tmpNewMolecule, tmpAdjList);
            List<IAtomContainer> tmpIsolatedRings = tmpRingSearch.isolatedRingFragments();
            List<IAtomContainer> tmpFusedRings = tmpRingSearch.fusedRingFragments();
            for(IAtomContainer tmpReferenceRing : this.ringSugars){
                for(IAtomContainer tmpIsolatedRing : tmpIsolatedRings){
                    if (this.univIsomorphismTester.isIsomorph(tmpReferenceRing, tmpIsolatedRing)) {
                        //do not remove rings with non-single exocyclic bonds
                        boolean tmpAreAllExocyclicBondsSingle = this.areAllExocyclicBondsSingle(tmpIsolatedRing, tmpNewMolecule);
                        if (!tmpAreAllExocyclicBondsSingle) {
                            continue;
                        }
                        //do not remove rings with 'too few' attached oxygens
                        int tmpExocyclicOxygenCount = this.getAttachedOxygenAtomCount(tmpIsolatedRing, tmpNewMolecule);
                        int tmpAtomsInRing = tmpIsolatedRing.getAtomCount();
                        boolean tmpAreEnoughOxygensAttached = this.doesRingHaveEnoughOxygenAtomsAttached(tmpAtomsInRing, tmpExocyclicOxygenCount);
                        if (!tmpAreEnoughOxygensAttached) {
                            continue;
                        }
                        //TODO: Move to constant
                        tmpNewMolecule.setProperty("CONTAINS_RING_SUGAR", 1);
                        //remove found ring
                        for(IAtom tmpAtom : tmpIsolatedRing.atoms()){
                            if (tmpNewMolecule.contains(tmpAtom)) {
                                tmpNewMolecule.removeAtom(tmpAtom);
                            }
                        }
                        tmpIsolatedRings.remove(tmpIsolatedRing);
                    }
                }
            }
            //***removing linear sugars***
            boolean tmpContainsLinearSugar = false;
            //iterating over linear sugar patterns
            for(DfPattern tmpPattern : this.linearSugarPatterns) {
                int[][] tmpUniqueMappings = tmpPattern.matchAll(tmpNewMolecule).uniqueAtoms().toArray();
                //Linear sugars that are detected as part of a larger ring should not be removed
                //iterating over mappings of one pattern
                for (int[] tmpMapping : tmpUniqueMappings) {
                    boolean tmpIsPartOfRing = false;
                    //iterating over atoms of one mapping
                    for (int i = 0; i < tmpMapping.length; i++) {
                        //TODO: Why this if?
                        if(tmpNewMolecule.getAtomCount() >= i) {
                            //iterating over fused rings to check whether atom is part of a ring
                            for(IAtomContainer tmpFusedRing : tmpFusedRings){
                                if(tmpFusedRing.contains(tmpNewMolecule.getAtom(i))) {
                                    //mapped atom is in a ring
                                    tmpIsPartOfRing = true;
                                    //breaking iteration of rings
                                    break;
                                }
                            }
                            if (tmpIsPartOfRing) {
                                //breaking iteration of atoms in mapping
                                break;
                            }
                            //doing the same for isolated rings
                            for(IAtomContainer tmpIsolatedRing : tmpIsolatedRings) {
                                if(tmpIsolatedRing.contains(tmpNewMolecule.getAtom(i))) {
                                    //mapped atom is in a ring
                                    tmpIsPartOfRing = true;
                                    break;
                                }
                            }
                            if (tmpIsPartOfRing) {
                                break;
                            }
                        }
                    }
                    if(!tmpIsPartOfRing ) {
                        ArrayList<IAtom> tmpAtomsToRemove = new ArrayList<>(tmpMapping.length);
                        for (int i = 0; i < tmpMapping.length; i++) {
                            //TODO: Why this if?
                            if(tmpNewMolecule.getAtomCount() >= i) {
                                try {
                                    IAtom atomToRemove = tmpNewMolecule.getAtom(i);
                                    tmpAtomsToRemove.add(atomToRemove);
                                } catch (IndexOutOfBoundsException anException) {
                                    SugarRemovalUtility.LOGGER.log(Level.SEVERE, anException.toString(), anException);
                                }
                            }
                        }
                        //TODO: To prevent removal of sugars that are not there anymore? Move to constant or option
                        if(tmpAtomsToRemove.size() >= 4) {
                            tmpContainsLinearSugar = true;
                        }
                        for(IAtom tmpAtom : tmpAtomsToRemove) {
                            if (tmpNewMolecule.contains(tmpAtom)) {
                                tmpNewMolecule.removeAtom(tmpAtom);
                            }
                        }
                    }
                }
            }
            //TODO: Move to constant
            if(tmpContainsLinearSugar) {
                tmpNewMolecule.setProperty("CONTAINS_LINEAR_SUGAR", 1);
            }
            //select only the biggest remaining unconnected part of the molecule
            tmpNewMolecule = SugarRemovalUtility.selectBiggestUnconnectedFragment(tmpNewMolecule);
        //TODO: maybe catch some exceptions in between already
        } catch (CloneNotSupportedException | ConcurrentModificationException | IndexOutOfBoundsException | CDKException anException) {
            SugarRemovalUtility.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            return null;
        }
        return tmpNewMolecule;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * TODO: Add doc, note that this method does not clear too small structures and may return null!
     *
     * @param aMolecule
     * @return
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
            //if something went wrong
            return null;
        }
        tmpBiggestFragment.setProperties(tmpProperties);
        return tmpBiggestFragment;
    }

    /**
     * TODO: Add doc, note that this method does not clear too small structures and may return null!
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
     * TODO
     */
    public static List<IAtomContainer> partitionAndSortUnconnectedFragments(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException {
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
     * TODO
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
     * order.
     * <p>
     *     The method iterates over all cyclic atoms and all of their bonds. So the runtime scales linear with the number
     *     of cyclic atoms and their connected bonds.
     * </p>
     *
     * @param aRingToTest the ring fragment to test; exocyclic bonds do not have to be included in the fragment but if it
     *                    is a fused system of multiple rings, the internal interconnecting bonds of the different rings
     *                    need to be included; all its atoms need to be exactly the same objects as in the second atom
     *                    container parameter
     * @param anOriginalMolecule the molecule that contains the ring under investigation; The exocyclic bonds will be
     *                           queried from it
     * @return true, if all exocyclic bonds connected to the ring are of single order
     * @throws NullPointerException if a parameter is 'null'
     */
    private boolean areAllExocyclicBondsSingle(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule) throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule is 'null'");
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

    //TODO/discuss: Include N-, S- and C-glycosidic bonds?
    //TODO/discuss: Include bonds that are not of type -X- but also of type -X(R)R etc., e.g. in adenosine?
    //Note: detects also ester and ether bonds which is not a bad thing because they occur frequently in NPs
    /**
     * TODO
     */
    private boolean hasGlycosidicBond(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule) throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule is 'null'");
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
                    boolean tmpIsOxygen = (tmpSymbol == "O");
                    if (tmpIsOxygen) {
                        List<IBond> tmpConnectedBondsList = anOriginalMolecule.getConnectedBondsList(tmpAtom);
                        boolean tmpHasOnlyTwoBonds = (tmpConnectedBondsList.size() == 2);
                        boolean tmpAllBondsAreSingle = (BondManipulator.getMaximumBondOrder(tmpConnectedBondsList) == IBond.Order.SINGLE);
                        if (tmpHasOnlyTwoBonds && tmpAllBondsAreSingle) {
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

    //TODO/discuss: Also count N and S (or all hetero atoms)?
    /**
     * Gives the number of attached exocyclic oxygen atoms of a given ring fragment of an original atom container.
     * <p>
     *     The method iterates over all cyclic atoms and all of their connected atoms. So the runtime scales linear with the number
     *     of cyclic atoms and their connected atoms.
     * </p>
     *
     * @param aRingToTest the ring fragment to test; exocyclic bonds do not have to be included in the fragment but if it
     *                    is a fused system of multiple rings, the internal interconnecting bonds of the different rings
     *                    need to be included; all its atoms need to be exactly the same objects as in the second atom
     *                    container parameter
     * @param anOriginalMolecule the molecule that contains the ring under investigation; The exocyclic bonds will be
     *                           queried from it
     * @return number of attached exocyclic oxygen atoms
     * @throws NullPointerException if a parameter is 'null'
     */
   private int getAttachedOxygenAtomCount(IAtomContainer aRingToTest, IAtomContainer anOriginalMolecule) throws NullPointerException {
        Objects.requireNonNull(aRingToTest, "Given ring atom container is 'null'");
        Objects.requireNonNull(anOriginalMolecule, "Given atom container representing the original molecule is 'null'");
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
                boolean tmpIsOxygen = tmpSymbol.matches("O");
                boolean tmpIsInRing = aRingToTest.contains(tmpConnectedAtom);
                if (tmpIsOxygen && !tmpIsInRing) {
                    tmpExocyclicOxygenCounter++;
                }
            }
        }
        return tmpExocyclicOxygenCounter;
    }

    //TODO: Revise doc and parameter tests
    /**
     * Simple decision making function for deciding whether a possible sugar ring has enough exocyclic oxygen atom
     * attached to it. This number should be higher or equal to the number of atoms in the ring (including the cyclic
     * oxygen atom) divided by two. So at least 3 attached exocyclic oxygen atoms for a six-membered ring, 2 for a
     * five-membered ring etc. A simple integer division is used.
     *
     * @param aNumberOfAtomsInRing number of atoms in the possible sugar ring, including the cyclic oxygen atom
     * @param aNumberOfAttachedExocyclicOxygenAtoms number of attached exocyclic oxygen atom of the ring under
     *                                              investigation
     * @return true, if the number of attached exocyclic oxygen atoms is at least half of the number of atoms in the ring
     */
    private boolean doesRingHaveEnoughOxygenAtomsAttached(int aNumberOfAtomsInRing, int aNumberOfAttachedExocyclicOxygenAtoms) {
        double tmpAttachedOxygensToAtomsInRingRatio = ((double) aNumberOfAttachedExocyclicOxygenAtoms / (double) aNumberOfAtomsInRing);
        boolean tmpMeetsThreshold = (tmpAttachedOxygensToAtomsInRingRatio >= this.attachedOxygensToAtomsInRingRatioThreshold);
        return tmpMeetsThreshold;
    }

    /**
     * TODO
     */
    private IAtomContainer clearTooSmallStructures(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            return aMolecule;
        }
        IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
            boolean tmpIsTooSmall = this.isTooSmall(tmpComponent);
            if (tmpIsTooSmall) {
                for (IAtom tmpAtom : tmpComponent.atoms()) {
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
     * TODO
     */
    private boolean isTooSmall(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        boolean tmpIsTooSmall;
        if (this.structuresToKeepMode == StructuresToKeepMode.ALL) {
            tmpIsTooSmall = false;
        } else if (this.structuresToKeepMode == StructuresToKeepMode.HEAVY_ATOM_COUNT) {
            int tmpHeavyAtomCount = AtomContainerManipulator.getHeavyAtoms(aMolecule).size();
            tmpIsTooSmall = tmpHeavyAtomCount < this.mwOrHacThreshold;
        } else if (this.structuresToKeepMode == StructuresToKeepMode.MOLECULAR_WEIGHT) {
            double tmpMolWeight = AtomContainerManipulator.getMass(aMolecule, AtomContainerManipulator.MolWeight);
            tmpIsTooSmall = tmpMolWeight < this.mwOrHacThreshold;
        } else {
            throw new UnsupportedOperationException("Undefined StructuresToKeepMode setting!");
        }
        return tmpIsTooSmall;
    }

    /**
     * TODO
     */
    private boolean isTerminal(IAtomContainer aSubstructure, IAtomContainer aParentMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
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
     * TODO
     */
    private IAtomContainer removeSugarCandidates(IAtomContainer aMolecule, List<IAtomContainer> aCandidateList) throws NullPointerException, IllegalArgumentException {
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
                        SugarRemovalUtility.LOGGER.log(Level.WARNING, aCloneNotSupportedException.toString(), aCloneNotSupportedException);
                        throw new IllegalArgumentException("Could not clone one candidate and therefore not determine " +
                                "whether it is terminal or not.");
                    }
                    if (tmpIsTerminal) {
                        for (IAtom tmpAtom : tmpCandidate.atoms()) {
                            if (tmpNewMolecule.contains(tmpAtom)) {
                                tmpNewMolecule.removeAtom(tmpAtom);
                            }
                        }
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
        //To clear away unconnected, too small structures and generate valid valences, the post-processing method must be called
        return tmpNewMolecule;
    }

    /**
     * TODO
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
     * TODO
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
                        boolean tmpAreEnoughOxygensAttached = this.doesRingHaveEnoughOxygenAtomsAttached(tmpAtomsInRing, tmpExocyclicOxygenCount);
                        if (!tmpAreEnoughOxygensAttached) {
                            continue;
                        }
                    }
                    //if sugar ring has not been excluded yet, the molecule contains sugars, although they might not be terminal
                    tmpSugarCandidates.add(tmpIsolatedRing);
                }
            }
        }
        return tmpSugarCandidates;
    }

    /**
     * TODO
     */
    private List<IAtomContainer> getLinearSugarCandidates(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        IAtomContainer tmpNewMolecule = aMolecule;
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(aMolecule.getAtomCount() / 2);
        HashSet<Integer> tmpSugarCandidateAtomsSet = new HashSet<>(aMolecule.getAtomCount() + 2, 1);
        for(DfPattern tmpLinearSugarPattern : this.linearSugarPatterns) {
            /*unique in this case means that the same match cannot be in this collection multiple times but they can
            still overlap! Overlapping atoms are removed in the following lines.*/
            //TODO/discuss: Is there a better way to get non-overlapping matches?
            Mappings tmpMappings = tmpLinearSugarPattern.matchAll(tmpNewMolecule);
            Mappings tmpUniqueMappings = tmpMappings.uniqueAtoms();
            Iterable<IAtomContainer> tmpUniqueSubstructureMappings = tmpUniqueMappings.toSubstructures();
            for (IAtomContainer tmpMatchedStructure : tmpUniqueSubstructureMappings) {
                for (IAtom tmpAtom : tmpMatchedStructure.atoms()) {
                    int tmpAtomIndex = tmpAtom.getProperty(SugarRemovalUtility.INDEX_PROPERTY_KEY);
                    boolean tmpIsAtomAlreadyInCandidates = tmpSugarCandidateAtomsSet.contains(tmpAtomIndex);
                    if (tmpIsAtomAlreadyInCandidates) {
                        tmpMatchedStructure.removeAtom(tmpAtom);
                    } else {
                        tmpSugarCandidateAtomsSet.add(tmpAtomIndex);
                    }
                }
                if (!tmpMatchedStructure.isEmpty()) {
                    tmpSugarCandidates.add(tmpMatchedStructure);
                }
            }
        }
        if (!this.removeLinearSugarsInRing && !tmpSugarCandidates.isEmpty()) {
            int[][] tmpAdjList = GraphUtil.toAdjList(tmpNewMolecule);
            RingSearch tmpRingSearch = new RingSearch(tmpNewMolecule, tmpAdjList);
            for (int i = 0; i < tmpSugarCandidates.size(); i++) {
                IAtomContainer tmpCandidate = tmpSugarCandidates.get(i);
                //TODO/discuss: discard whole candidate substructure if at least one of its atoms is in a ring (as done now)
                // or remove only the atoms that are in a ring from the candidate substructure?
                boolean tmpIsInRing = false;
                for (IAtom tmpAtom : tmpCandidate.atoms()) {
                    if (tmpRingSearch.cyclic(tmpAtom)) {
                        tmpIsInRing = true;
                        break;
                    }
                }
                if (tmpIsInRing) {
                    tmpSugarCandidates.remove(i);
                    //The removal shifts the remaining indices!
                    i = i - 1;
                }
            }
        }
        return tmpSugarCandidates;
    }
    //</editor-fold>
}