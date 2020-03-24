/**
 * TODO: Add license header
 */
package de.unijena.cheminf.deglycosylation;

/**
 * TODO:
 * - If all remaining structures are kept, hydroxy groups of the sugars are also kept. If only terminal sugars should be
 * removed, no sugar is removed in this case! So remove sugars with attached hydroxy groups also?
 * - implement removal of linear sugars and of both circular and linear
 * - implement hasXy methods
 * - Note: The situation at the previous connection point is unclear. When circular sugars are removed, a hydroxy
 * group remains at the core structure (if there was a glycosidic bond). But for the linear sugars, no general statement
 * like this can be made.
 * - add hetero atom count as StructureToKeepMode option?
 * - add detection of glycosidic bond for linear sugars?
 * - see all the to dos
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
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;

import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * TODO: Add doc
 */
public final class SugarRemovalUtility {
    //<editor-fold desc="Enum StructuresToKeepMode">
    /**
     * TODO
     */
    public static enum StructuresToKeepMode {
        /**
         * TODO
         */
        ALL (0),

        /**
         * TODO
         */
        HEAVY_ATOM_COUNT (5),

        /**
         * TODO
         */
        MOLECULAR_WEIGHT (60);

        /**
         * TODO
         */
        private final int defaultThreshold;

        /**
         * TODO
         */
        StructuresToKeepMode(int aDefaultValue) {
            this.defaultThreshold = aDefaultValue;
        }

        /**
         * TODO
         */
        public int getDefaultThreshold() {
            return this.defaultThreshold;
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final constants">
    /**
     * TODO
     */
    public static final String CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY = "CONTAINS_CIRCULAR_SUGAR";

    /**
     * TODO
     */
    public static final String CONTAINS_LINEAR_SUGAR_PROPERTY_KEY = "CONTAINS_LINEAR_SUGAR";

    /**
     * TODO
     */
    public static final String INDEX_PROPERTY_KEY = "SUGAR_REMOVAL_UTILITY_INDEX";

    /**
     * TODO: Add doc, add names of sugars
     */
    public static final String[] LINEAR_SUGARS_SMILES = {
            "C(C(C(C(C(C=O)O)O)O)O)O",
            "C(C(CC(C(CO)O)O)O)(O)=O",
            "C(C(C(CC(=O)O)O)O)O",
            "C(C(C(C(C(CO)O)O)O)=O)O",
            "C(C(C(C(C(CO)O)O)O)O)O",
            "C(C(C(C(CC=O)O)O)O)O",
            "OCC(O)C(O)C(O)C(O)CO",
            "O=CC(O)C(O)C(O)C(O)CO",
            "CCCCC(O)C(=O)O", //TODO: Is this a sugar?
            "CC(=O)CC(=O)CCC(=O)O", //TODO: Is this a sugar?
            "O=C(O)CC(O)CC(=O)O", //TODO: Is this a sugar?
            "O=C(O)C(=O)C(=O)C(O)C(O)CO",
            "O=C(O)CCC(O)C(=O)O", //TODO: Is this a sugar?
            "O=CC(O)C(O)C(O)C(O)CO",
            "O=C(CO)C(O)C(O)CO"};

    /**
     * TODO: Add doc, add names/description
     */
    public static final String [] RING_SUGARS_SMILES = {
            "C1CCOC1",
            "C1CCOCC1",
            "C1CCCOCC1"};

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
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * TODO: Add doc, add parameters?
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
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set/add/clear">
    /**
     * TODO: Add doc, add test to check whether given new structure is isomorph to an already present structure
     */
    public void addCircularSugar(IAtomContainer aCircularSugar) throws NullPointerException, IllegalArgumentException, IllegalStateException {
        Objects.requireNonNull(aCircularSugar, "Given atom container is 'null'");
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
        Objects.requireNonNull(aLinearSugar, "Given atom container is 'null'");
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
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * TODO
     */
    public IAtomContainer removeCircularSugars(IAtomContainer aMolecule, boolean aShouldBeCloned) throws NullPointerException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        this.setIndices(tmpNewMolecule);
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
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (tmpContainsSugar) {
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY, true);
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
                        boolean tmpIsTerminal = this.isTerminal(tmpCandidate, tmpNewMolecule);
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
                            this.clearTooSmallStructures(tmpNewMolecule);
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
                for (IAtomContainer tmpSugarRing : tmpSugarCandidates) {
                    for (IAtom tmpAtom : tmpSugarRing.atoms()) {
                        if (tmpNewMolecule.contains(tmpAtom)) {
                            tmpNewMolecule.removeAtom(tmpAtom);
                        }
                    }
                }
            }
            //if too small / too light, unconnected structures should be discarded, this is done now
            //otherwise, the possibly unconnected atom container is returned
            //If only terminal sugars are removed, the resulting structure may still be too small to keep!
            if (this.structuresToKeepMode != StructuresToKeepMode.ALL) {
                this.clearTooSmallStructures(tmpNewMolecule);
            }
            //post-processing
            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpNewMolecule);
                CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(tmpNewMolecule);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
            }
        }
        //May be empty and may be unconnected, based on the settings
        return tmpNewMolecule;
    }

    /**
     * TODO
     */
    public IAtomContainer removeLinearSugars(IAtomContainer aMolecule, boolean aShouldBeCloned) throws NullPointerException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        IAtomContainer tmpNewMolecule;
        if (aShouldBeCloned) {
            tmpNewMolecule = aMolecule.clone();
        } else {
            tmpNewMolecule = aMolecule;
        }
        this.setIndices(tmpNewMolecule);
        List<IAtomContainer> tmpSugarCandidates = new ArrayList<>(aMolecule.getAtomCount() / 2);
        HashSet<Integer> tmpSugarCandidateAtomsSet = new HashSet<>(aMolecule.getAtomCount() + 2, 1);
        for(DfPattern tmpLinearSugarPattern : this.linearSugarPatterns) {
            /*unique in this case means that the same match cannot be in this collection multiple times but they can
            still overlap! Overlapping atoms are removed in the following lines.*/
            //TODO: Is there a better way to get non-overlapping matches?
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
        //note: this does only mean that there are matches of the linear sugar patterns! They might not be terminal or
        //might be contained in a ring!
        boolean tmpContainsSugar = !tmpSugarCandidates.isEmpty();
        if (tmpContainsSugar) {
            tmpNewMolecule.setProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY, true);
            if (!this.removeLinearSugarsInRing) {
                int[][] tmpAdjList = GraphUtil.toAdjList(tmpNewMolecule);
                RingSearch tmpRingSearch = new RingSearch(tmpNewMolecule, tmpAdjList);
                for (int i = 0; i < tmpSugarCandidates.size(); i++) {
                    IAtomContainer tmpCandidate = tmpSugarCandidates.get(i);
                    //TODO: discard whole candidate substructure if at least one of its atoms is in a ring (as done now)
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
                if (tmpSugarCandidates.isEmpty()) {
                    //to skip the remaining code and exit the method
                    tmpContainsSugar = false;
                }
            }
        }
        if (tmpContainsSugar) {
            //TODO: duplicate!
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
                        boolean tmpIsTerminal = this.isTerminal(tmpCandidate, tmpNewMolecule);
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
                            this.clearTooSmallStructures(tmpNewMolecule);
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
            //TODO: Post-processing is a duplicate, could be moved to its own method
            //post-processing
            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpNewMolecule);
                CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(tmpNewMolecule);
            } catch (CDKException aCDKException) {
                SugarRemovalUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
            }
        }
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
    //TODO: Add method to select heaviest, unconnected structure
    /**
     * TODO: Add doc
     *
     * @param aMolecule
     * @return
     */
    public static IAtomContainer selectBiggestUnconnectedFragment(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        Map<Object, Object> tmpProperties = aMolecule.getProperties();
        IAtomContainerSet tmpUnconnectedFragments = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        if(tmpUnconnectedFragments != null && tmpUnconnectedFragments.getAtomContainerCount() >= 1) {
            IAtomContainer tmpBiggestFragment = tmpUnconnectedFragments.getAtomContainer(0);
            for(IAtomContainer tmpFragment : tmpUnconnectedFragments.atomContainers()){
                if(tmpFragment.getAtomCount() > tmpBiggestFragment.getAtomCount()){
                    tmpBiggestFragment = tmpFragment;
                }
            }
            aMolecule = tmpBiggestFragment;
            int tmpHeavyAtomCount = 0;
            for(IAtom tmpAtom : aMolecule.atoms()){
                if(!tmpAtom.getSymbol().equals("H")){
                    tmpHeavyAtomCount++;
                }
            }
            //Discard 'too small' remaining fragments
            //TODO: Move to constants or rather option or remove!
            if(tmpHeavyAtomCount < 5){
                return null;
            }
        } else {
            return null;
        }
        aMolecule.setProperties(tmpProperties);
        return aMolecule;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * TODO
     */
    private void setIndices(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int tmpIndex = 0;
        for (IAtom tmpAtom : aMolecule.atoms()) {
            tmpAtom.setProperty(SugarRemovalUtility.INDEX_PROPERTY_KEY, tmpIndex);
            tmpIndex++;
        }
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
            //TODO: 'NoSuchAtomException: Atom is not a member of this AtomContainer' is thrown sometimes in the following line
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

    //TODO: Include N-, S- and C-glycosidic bonds?
    // TODO: Include bonds that are not of type -X- but also of type -X(R)R etc., e.g. in adenosine?
    //Note: detects also ester bonds which is not a bad thing because they occur frequently in NPs
    /**
     * TODO!
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

    //TODO: Also count N and S (or all hetero atoms)?
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
     *TODO
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
        //TODO: Is there a better way?
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
    //</editor-fold>
}
