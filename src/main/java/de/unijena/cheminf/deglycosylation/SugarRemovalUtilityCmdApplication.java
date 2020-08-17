/*
 * MIT License
 *
 * Copyright (c) 2020 Jonas Schaub, Achim Zielesny, Christoph Steinbeck, Maria Sorokina
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.deglycosylation;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.cdkbook.SMILESFormatMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.IChemFormatMatcher;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Objects;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * Controller of the Sugar Removal Utility command line application. It can be used to remove sugar moieties from
 * molecules in a given data set, according to "Schaub J, Zielesny A, Steinbeck C, Sorokina M.
 * Too sweet: cheminformatics for deglycosylation in natural products, 02 August 2020, PREPRINT (Version 1)", available
 * at Research Square [<a href="https://doi.org/10.21203/rs.3.rs-50194/v1">DOI:10.21203/rs.3.rs-50194/v1</a>]. This class
 * basically instantiates the SugarRemovalUtility class with the settings specified in the command line arguments and
 * uses it to iterate over all molecules found in the given file and remove their sugar moieties. Also, a file detailing
 * the deglycosylated cores and removed sugar moieties for each molecule is written as output.
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class SugarRemovalUtilityCmdApplication {
    //<editor-fold desc="Private variables">
    /**
     * SugarRemovalUtility instance used to detect and remove the sugar moieties.
     */
    private SugarRemovalUtility sugarRemovalUtil;

    /**
     * Setting that specifies whether to remove circular, linear, or both types of moieties from the molecules.
     */
    private int typeOfMoietiesToRemove;

    /**
     * Input file containing molecule set, either MDL Molfile, MDL Structure data file (SDF), or SMILES file
     */
    private File inputFile;
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Character used to separate the values in the output CSV file.
     */
    private static final String OUTPUT_FILE_SEPARATOR = ";";

    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(SugarRemovalUtilityCmdApplication.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * Constructor that only needs the input file and the type of moieties to remove. All SugarRemovalUtility settings
     * will be set to their default values. Instantiates this class and the SugarRemovalUtility.
     *
     * @param aFile input file containing a set of molecules, either MDL Molfile, MDL Structure data file (SDF), or
     *              SMILES file (of format: [SMILES string][space][name] in each line, see example file); the given file
     *              must exist and be accessible and readable; the file type extension is not important for the
     *              determination of the file type; the final test for whether the file is suitable is done in execute()
     * @param aTypeOfMoietiesToRemove an int of value ["1","2","3"] indicating which type of sugar moieties should be removed,
     *                                "1" for circular sugar moieties, "2" for linear sugar moieties, or "3" for
     *                                circular AND linear sugar moieties; check isLegalTypeOfMoietiesToRemove(int) for
     *                                the correct mapping of value to option
     * @throws IllegalArgumentException if any parameter does not meet the specified requirements
     */
    public SugarRemovalUtilityCmdApplication(File aFile, int aTypeOfMoietiesToRemove) throws IllegalArgumentException {
        this(
                aFile,
                aTypeOfMoietiesToRemove,
                SugarRemovalUtility.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_DEFAULT,
                SugarRemovalUtility.REMOVE_ONLY_TERMINAL_SUGARS_DEFAULT,
                SugarRemovalUtility.PRESERVATION_MODE_DEFAULT.ordinal(),
                SugarRemovalUtility.PRESERVATION_MODE_DEFAULT.getDefaultThreshold(),
                SugarRemovalUtility.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_DEFAULT,
                SugarRemovalUtility.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT,
                SugarRemovalUtility.DETECT_LINEAR_SUGARS_IN_RINGS_DEFAULT,
                SugarRemovalUtility.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_DEFAULT,
                SugarRemovalUtility.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_DEFAULT,
                SugarRemovalUtility.DETECT_LINEAR_ACIDIC_SUGARS_DEFAULT,
                SugarRemovalUtility.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_DEFAULT
        );
    }

    /**
     * Constructor that needs all settings explicitly specified. Instantiates this class and the SugarRemovalUtility with
     * the given settings.
     *
     * @param aFile input file containing a set of molecules, either MDL Molfile, MDL Structure data file (SDF), or
     *              SMILES file (of format: [SMILES string][space][name] in each line, see example file); the given file
     *              must exist and be accessible and readable; the file type extension is not important for the
     *              determination of the file type; the final test for whether the file is suitable is done in execute()
     * @param aTypeOfMoietiesToRemove an int of value ["1","2","3"] indicating which type of sugar moieties should be removed,
     *                                "1" for circular sugar moieties, "2" for linear sugar moieties, or "3" for
     *                                circular AND linear sugar moieties; check isLegalTypeOfMoietiesToRemove(int) for
     *                                the correct mapping of value to option
     * @param aDetectCircularSugarsOnlyWithOGlycosidicBondSetting true if circular sugars should be detected (and
     *                                                            removed) only if they have an O-glycosidic bond to
     *                                                            another moiety or the core of the molecule
     * @param aRemoveOnlyTerminalSugarsSetting true if only terminal sugar moieties should be removed; if this setting
     *                                         is set to "true", the input molecules must all consist of one connected
     *                                         structure, respectively; if they already contain multiple, disconnected
     *                                         structures (e.g. counter-ions), the respective molecules are ignored.
     * @param aStructureToKeepModeSetting an int of value ["0","1","2"] indicating which preservation mode to use; This
     *                                    specifies under what circumstances to discard structures that get disconnected
     *                                    from the central core in the sugar removal process, "0" to preserve all
     *                                    disconnected structures (note: this might lead to no circular sugar moieties
     *                                    being detected, depending on the other settings), "1" to remove disconnected
     *                                    structures that do not have enough heavy atoms, or "2" to remove disconnected
     *                                    structures that do not have a sufficient molecular weight; check
     *                                    SugarRemovalUtility.PreservationModeOption enum for the correct mapping of
     *                                    value to option, it corresponds to the ordinal values of the enum constants
     * @param aStructureToKeepModeThresholdSetting an int specifying the threshold of the preservation mode, i.e. how
     *                                             many heavy atoms a disconnected structure needs to have at least to
     *                                             be not removed or how heavy (in terms of its molecular weight) it
     *                                             needs to be; the value must be positive; if aStructureToKeepModeSetting
     *                                             was passed the value "0" (preserve all structures), this argument must
     *                                             also be passed a zero value; in the opposite case, this argument must
     *                                             be passed a non-zero value if aStructureToKeepModeSetting was given
     *                                             the value 1 or 2
     * @param aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting true if circular sugars should be detected
     *                                                                       (and removed) only if they have a sufficient
     *                                                                       number of attached exocyclic oxygen atoms
     * @param anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting a double value giving the minimum attached
     *                                                                 exocyclic oxygen atoms to atom number in the ring
     *                                                                 ratio a circular sugar needs to have to be
     *                                                                 detected as such, e.g 0.5 where a six-membered
     *                                                                 ring needs to have at least 3 oxygen atoms
     *                                                                 attached; if
     *                                                                 aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
     *                                                                 was passed the value "false" (detect circular sugars
     *                                                                 neglecting their number of attached exocyclic oxygen
     *                                                                 atoms), this argument must be passed a zero
     *                                                                 value; in the opposite case, this argument must
     *                                                                 be passed a non-zero value if the other
     *                                                                 argument was given the value "true"; the number
     *                                                                 must be positive.
     * @param aDetectLinearSugarsInRingsSetting true if linear sugars in rings should be detected (and removed)
     * @param aLinearSugarCandidateMinSizeSetting an int specifying the minimum number of carbon atoms a linear sugar
     *                                            needs to have to be detected as such; the value must be positive and
     *                                            higher than or equal to 1 and also smaller than the linear sugar
     *                                            candidate maximum size
     * @param aLinearSugarCandidateMaxSizeSetting an int specifying the maximum number of carbon atoms a linear sugar
     *                                            needs to have to be detected as such; the integer number must be
     *                                            positive and higher than or equal to 1 and also higher than the linear
     *                                            sugar candidate minimum size
     * @param aDetectLinearAcidicSugarsSetting true if linear acidic sugars should be included in the set of linear
     *                                         sugar patterns for the initial detection
     * @param aDetectSpiroRingsAsCircularSugarsSetting true if spiro rings (rings that share one atom with another cycle)
     *                                                 should be included in the circular sugar detection
     * @throws IllegalArgumentException if any parameter does not meet the specified requirements
     */
    public SugarRemovalUtilityCmdApplication(File aFile,
                                             int aTypeOfMoietiesToRemove,
                                             boolean aDetectCircularSugarsOnlyWithOGlycosidicBondSetting,
                                             boolean aRemoveOnlyTerminalSugarsSetting,
                                             int aStructureToKeepModeSetting,
                                             int aStructureToKeepModeThresholdSetting,
                                             boolean aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting,
                                             double anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting,
                                             boolean aDetectLinearSugarsInRingsSetting,
                                             int aLinearSugarCandidateMinSizeSetting,
                                             int aLinearSugarCandidateMaxSizeSetting,
                                             boolean aDetectLinearAcidicSugarsSetting,
                                             boolean aDetectSpiroRingsAsCircularSugarsSetting)
            throws IllegalArgumentException {
        if (!aFile.exists() || !aFile.canRead() || !aFile.isFile()) {
            throw new IllegalArgumentException("Given path to SMILES, SD, or MOL file does not exist or " +
                    "the file cannot be read.");
        }
        this.inputFile = aFile;
        //determination of file type is done in execute(), might cause exceptions there!
        if (!SugarRemovalUtilityCmdApplication.isLegalTypeOfMoietiesToRemove(aTypeOfMoietiesToRemove)) {
            throw new IllegalArgumentException("Type of moieties to remove must be either 1 (circular sugar moieties), " +
                    "2 (linear sugar moieties), or 3 (both).");
        }
        this.typeOfMoietiesToRemove = aTypeOfMoietiesToRemove;
        this.sugarRemovalUtil = new SugarRemovalUtility();
        this.sugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(aDetectCircularSugarsOnlyWithOGlycosidicBondSetting);
        this.sugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(aRemoveOnlyTerminalSugarsSetting);
        SugarRemovalUtility.PreservationModeOption tmpStructureToKeepMode = null;
        for (SugarRemovalUtility.PreservationModeOption tmpEnumConstant : SugarRemovalUtility.PreservationModeOption.values()) {
            if (aStructureToKeepModeSetting == tmpEnumConstant.ordinal()) {
                tmpStructureToKeepMode = tmpEnumConstant;
            }
        }
        if (Objects.isNull(tmpStructureToKeepMode)) {
            throw new IllegalArgumentException("Given structure to keep mode setting does not correspond to an ordinal " +
                    "in the StructureToKeepModeOption enumeration.");
        }
        this.sugarRemovalUtil.setPreservationModeSetting(tmpStructureToKeepMode);
        if (aStructureToKeepModeThresholdSetting < 0) {
            throw new IllegalArgumentException("The structure to keep mode setting threshold " +
                    "cannot be negative.");
        }
        //case 1: Structure to keep mode setting is 'all' (ordinal 0). Then, a 0 threshold should be passed.
        // case 2: Structure to keep mode is not 'all' (ordinal nonzero). Then, a nonzero threshold should be passed.
        if (tmpStructureToKeepMode == SugarRemovalUtility.PreservationModeOption.ALL && aStructureToKeepModeThresholdSetting != 0) {
            throw new IllegalArgumentException("Structure to keep mode 'all' was selected. Therefore, passing a nonzero " +
                    "threshold makes no sense.");
        } else if (tmpStructureToKeepMode != SugarRemovalUtility.PreservationModeOption.ALL && aStructureToKeepModeThresholdSetting == 0) {
            throw new IllegalArgumentException("Please select structure to keep mode 'all' if the associated " +
                    "threshold should be 0.");
        }
        this.sugarRemovalUtil.setPreservationModeThresholdSetting(aStructureToKeepModeThresholdSetting);
        this.sugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting);
        boolean tmpIsFinite = Double.isFinite(anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting);
        boolean tmpIsNegative = (anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting < 0);
        if (!tmpIsFinite || tmpIsNegative) {
            throw new IllegalArgumentException("Given exocyclic oxygen atoms to atoms in ring ration threshold is NaN, " +
                    "infinite or negative.");
        }
        //case 1: If the number of exocyclic oxygen atoms is neglected, a 0 threshold should be passed
        // case 2: If it is detected, a nonzero threshold should be passed
        if (!aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
                && anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting != 0) {
            throw new IllegalArgumentException("The number of exocyclic oxygen atoms of circular sugars is neglected at detection. " +
                    "Therefore, passing a nonzero threshold makes no sense.");
        } else if (aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
                && anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting == 0) {
            throw new IllegalArgumentException("Please select to neglect the number of exocyclic oxygen atoms if the associated " +
                    "threshold should be 0.");
        }
        this.sugarRemovalUtil.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(anExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting);
        this.sugarRemovalUtil.setDetectLinearSugarsInRingsSetting(aDetectLinearSugarsInRingsSetting);
        if (aLinearSugarCandidateMinSizeSetting < 1 || aLinearSugarCandidateMaxSizeSetting < 1) {
            throw new IllegalArgumentException("The linear sugar candidate minimum and maximum sizes must both be equal or " +
                    "higher than 1.");
        }
        if (aLinearSugarCandidateMinSizeSetting > aLinearSugarCandidateMaxSizeSetting) {
            throw new IllegalArgumentException("The linear sugar candidate minimum size must be smaller " +
                    "than the maximum size.");
        }
        this.sugarRemovalUtil.setLinearSugarCandidateMinSizeSetting(aLinearSugarCandidateMinSizeSetting);
        this.sugarRemovalUtil.setLinearSugarCandidateMaxSizeSetting(aLinearSugarCandidateMaxSizeSetting);
        this.sugarRemovalUtil.setDetectLinearAcidicSugarsSetting(aDetectLinearAcidicSugarsSetting);
        this.sugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(aDetectSpiroRingsAsCircularSugarsSetting);
        this.sugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Executes the application. First, opens the given input file and determines its type. Then creates an output file
     * and a log file in the same directory as the input file. All settings are printed to console. The molecule set is
     * then iterated and all detected sugar moieties removed according to the given settings. All results are printed in
     * the output file and a small summary is printed to console at the end.
     *
     * @throws IOException if there are problems with reading and writing files
     * @throws SecurityException if the application is forbidden to access required files
     * @throws IllegalArgumentException if the given input file cannot be used or another setting in this class is
     * erroneous
     */
    public void execute() throws IOException, SecurityException, IllegalArgumentException {
        System.out.println("Trying to load file at " + this.inputFile.getAbsolutePath());
        FileInputStream tmpFileInputStream = new FileInputStream(this.inputFile);
        InputStreamReader tmpInputStreamReader = new InputStreamReader(tmpFileInputStream);
        BufferedReader tmpBuffReader = new BufferedReader(tmpInputStreamReader);
        FileReader tmpFileReader = new FileReader(this.inputFile);
        FormatFactory tmpCDKFormatFactory = new FormatFactory();
        tmpCDKFormatFactory.registerFormat((IChemFormatMatcher) SMILESFormatMatcher.getInstance());
        IChemFormat tmpFormat = tmpCDKFormatFactory.guessFormat(tmpBuffReader);
        if (Objects.isNull(tmpFormat)) {
            throw new IllegalArgumentException("Given file format cannot be determined.");
        }
        String tmpFormatClassName = tmpFormat.getFormatName();
        IIteratingChemObjectReader<IAtomContainer> tmpReader = null;
        switch (tmpFormatClassName) {
            case "MDL Structure-data file":
            case "MDL Molfile":
            case "MDL Molfile V2000":
            case "MDL Mol/SDF V3000":
                tmpReader = new IteratingSDFReader(tmpFileReader, DefaultChemObjectBuilder.getInstance(), true);
                System.out.println("Found MDL molfile / structure data file.");
                break;
            case "SMILES":
                tmpReader = new IteratingSMILESReader(tmpFileReader, DefaultChemObjectBuilder.getInstance());
                System.out.println("Found SMILES file.");
                break;
            default:
                throw new IllegalArgumentException("Given file type " + tmpFormatClassName + " cannot be used");
        }
        String tmpOutputRootDirectoryPath = this.inputFile.getAbsoluteFile().getParent() + File.separator;
        String tmpOutputFilePath = tmpOutputRootDirectoryPath + "deglycosylation_results.txt";
        File tmpOutputFile = new File(tmpOutputFilePath);
        FileWriter tmpOutputFileWriter = new FileWriter(tmpOutputFile);
        PrintWriter tmpOutputFilePrintWriter = new PrintWriter(tmpOutputFileWriter, true);
        String tmpOutputFileHeader = "Nr" + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR
                + "ID" + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR
                + "originalMoleculeSMILES" + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR
                + "deglycosylatedMoleculeSMILES" + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR
                + "hadOrHasSugars" + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR
                + "SugarMoietySMILES";
        tmpOutputFilePrintWriter.println(tmpOutputFileHeader);
        tmpOutputFilePrintWriter.flush();
        System.out.println("Results will be written to: " + tmpOutputFile.getAbsolutePath());
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputRootDirectoryPath + "Log.txt");
        } catch (IOException anIOException) {
            SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger tmpRootLogger = Logger.getLogger("");
        /*for (Handler tmpHandler : tmpRootLogger.getHandlers()) {
            tmpRootLogger.removeHandler(tmpHandler);
        }*/
        tmpRootLogger.addHandler(tmpLogFileHandler);
        tmpRootLogger.setLevel(Level.WARNING);
        System.out.println("Logging will be done in file at: " + tmpOutputRootDirectoryPath + "Log.txt");
        System.out.println("Given settings for Sugar Removal Utility:");
        switch (this.typeOfMoietiesToRemove) {
            case 1:
                System.out.println("Type of moieties to remove: circular sugars (1)");
                break;
            case 2:
                System.out.println("Type of moieties to remove: linear sugars (2)");
                break;
            case 3:
                System.out.println("Type of moieties to remove: circular AND linear sugars (3)");
                break;
            default:
                throw new IllegalArgumentException("Type of moieties to remove must be either 1 (circular sugar moieties), " +
                        "2 (linear sugar moieties), or 3 (both).");
        }
        System.out.println("Detect circular sugars only with O-glycosidic bond: "
                + this.sugarRemovalUtil.areOnlyCircularSugarsWithOGlycosidicBondDetected());
        System.out.println("Remove only terminal sugar moieties: " + this.sugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println("Structure to keep mode setting: " + this.sugarRemovalUtil.getPreservationModeSetting().name()
                + " (" + this.sugarRemovalUtil.getPreservationModeSetting().ordinal() + ")");
        System.out.println("Structure to keep mode threshold: "
                + this.sugarRemovalUtil.getPreservationModeThresholdSetting());
        System.out.println("Detect circular sugars only with enough exocyclic oxygen atoms: "
                + this.sugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected());
        System.out.println("Exocyclic oxygen atoms to atoms in ring ratio threshold: "
                + this.sugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        System.out.println("Detect linear sugars in rings: " + this.sugarRemovalUtil.areLinearSugarsInRingsDetected());
        System.out.println("Linear sugar candidate minimum size: "
                + this.sugarRemovalUtil.getLinearSugarCandidateMinSizeSetting() + " carbon atoms");
        System.out.println("Linear sugar candidate maximum size: "
                + this.sugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting() + " carbon atoms");
        System.out.println("Detect linear acidic sugars: " + this.sugarRemovalUtil.areLinearAcidicSugarsDetected());
        System.out.println("Detect spiro rings as circular sugars: "
                + this.sugarRemovalUtil.areSpiroRingsDetectedAsCircularSugars());
        System.out.println("Entering iteration of molecules...");
        this.sugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        int tmpMoleculeCounter = 1;
        int tmpFatalExceptionCounter = 0;
        int tmpMinorExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        int tmpContainsLinearSugarsCounter = 0;
        int tmpContainsCircularSugarsCounter = 0;
        int tmpBasicallyASugarCounter = 0;
        IAtomContainer tmpMolecule;
        String tmpID;
        String tmpOriginalMoleculeSMILES;
        boolean tmpHadSugars;
        String tmpDeglycosylatedMoleculeSMILES;
        List<IAtomContainer> tmpSugarMoieties;
        String tmpOutput = "";
        while (tmpReader.hasNext()) {
            try {
                tmpOutput = tmpMoleculeCounter + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR;
                try {
                    tmpMolecule = tmpReader.next();
                } catch (Exception anException) {
                    SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE,
                            anException.toString() + " Molecule number: " + tmpMoleculeCounter,
                            anException);
                    tmpLogFileHandler.flush();
                    tmpMolecule = null;
                }
                if (Objects.isNull(tmpMolecule) || tmpMolecule.isEmpty()) {
                    tmpOutput = tmpOutput.concat("[Molecule could not be read]");
                    tmpOutputFilePrintWriter.println(tmpOutput);
                    tmpOutputFilePrintWriter.flush();
                    tmpFatalExceptionCounter++;
                    continue;
                }
                if (!Objects.isNull(tmpMolecule.getProperty("Name"))) {
                    tmpID = tmpMolecule.getProperty("Name");
                } else if (!Objects.isNull(tmpMolecule.getProperty("name"))) {
                    tmpID = tmpMolecule.getProperty("name");
                } else if (!Objects.isNull(tmpMolecule.getProperty(CDKConstants.TITLE))) {
                    tmpID = tmpMolecule.getProperty(CDKConstants.TITLE);
                } else if (!Objects.isNull(tmpMolecule.getTitle())) {
                    tmpID = tmpMolecule.getTitle();
                } else if (!Objects.isNull(tmpMolecule.getID())) {
                    tmpID = tmpMolecule.getID();
                } else {
                    tmpID = "[No title available]";
                }
                if (tmpID.isEmpty() || tmpID.isBlank()) {
                    tmpID = "[No ID available]";
                }
                tmpOutput = tmpOutput.concat(tmpID + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR);
                try {
                    tmpOriginalMoleculeSMILES = tmpSmiGen.create(tmpMolecule);
                } catch (CDKException aCDKException) {
                    SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE,
                            aCDKException.toString() + " Molecule id: " + tmpID,
                            aCDKException);
                    tmpLogFileHandler.flush();
                    tmpOriginalMoleculeSMILES = "[SMILES code could not be generated]";
                    tmpMinorExceptionsCounter++;
                }
                tmpOutput = tmpOutput.concat(tmpOriginalMoleculeSMILES + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR);
                if (this.sugarRemovalUtil.areOnlyTerminalSugarsRemoved()) {
                    boolean tmpIsConnected = ConnectivityChecker.isConnected(tmpMolecule);
                    if (!tmpIsConnected) {
                        tmpDeglycosylatedMoleculeSMILES = "[Deglycosylation not possible because structure is unconnected]";
                        tmpOutput = tmpOutput.concat(tmpDeglycosylatedMoleculeSMILES);
                        tmpOutputFilePrintWriter.println(tmpOutput);
                        tmpOutputFilePrintWriter.flush();
                        tmpFatalExceptionCounter++;
                        continue;
                    }
                }
                try {
                    switch (this.typeOfMoietiesToRemove) {
                        case 1:
                            tmpSugarMoieties = this.sugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, false);
                            break;
                        case 2:
                            tmpSugarMoieties = this.sugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, false);
                            break;
                        case 3:
                            tmpSugarMoieties = this.sugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpMolecule, false);
                            break;
                        default:
                            throw new IllegalArgumentException("Type of moieties to remove must be either 1 (circular sugar moieties), " +
                                    "2 (linear sugar moieties), or 3 (both).");
                    }
                } catch (CloneNotSupportedException aCloneNotSupportedException) {
                    SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE,
                            aCloneNotSupportedException.toString() + " Molecule id: " + tmpID,
                            aCloneNotSupportedException);
                    tmpLogFileHandler.flush();
                    tmpDeglycosylatedMoleculeSMILES = "[Deglycosylation not possible because structure could not be cloned]";
                    tmpOutput = tmpOutput.concat(tmpDeglycosylatedMoleculeSMILES);
                    tmpOutputFilePrintWriter.println(tmpOutput);
                    tmpOutputFilePrintWriter.flush();
                    tmpFatalExceptionCounter++;
                    continue;
                }
                if (Objects.isNull(tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY))) {
                    tmpHadSugars = false;
                } else {
                    tmpHadSugars = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY);
                }
                if (tmpHadSugars) {
                    tmpSugarContainingMoleculesCounter++;
                }
                if (!Objects.isNull(tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY))) {
                    if (tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY)) {
                        tmpContainsLinearSugarsCounter++;
                    }
                }
                if (!Objects.isNull(tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY))) {
                    if (tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY)) {
                        tmpContainsCircularSugarsCounter++;
                    }
                }
                IAtomContainer tmpDeglycosylatedCore = tmpSugarMoieties.remove(0);
                if (tmpDeglycosylatedCore.isEmpty()) {
                    tmpDeglycosylatedMoleculeSMILES = "[empty]";
                    tmpBasicallyASugarCounter++;
                } else {
                    try {
                        tmpDeglycosylatedMoleculeSMILES = tmpSmiGen.create(tmpDeglycosylatedCore);
                    } catch (CDKException aCDKException) {
                        SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE,
                                aCDKException.toString() + " Molecule id: " + tmpID,
                                aCDKException);
                        tmpLogFileHandler.flush();
                        tmpDeglycosylatedMoleculeSMILES = "[SMILES code could not be generated]";
                        tmpMinorExceptionsCounter++;
                    }
                }
                tmpOutput = tmpOutput.concat(tmpDeglycosylatedMoleculeSMILES
                        + SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR
                        + tmpHadSugars);
                if (!tmpSugarMoieties.isEmpty()) {
                    for (IAtomContainer tmpMoiety : tmpSugarMoieties) {
                        String tmpSMILEScode = null;
                        try {
                            tmpSMILEScode = tmpSmiGen.create(tmpMoiety);
                        } catch (CDKException aCDKException) {
                            SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE,
                                    aCDKException.toString() + " Molecule id: " + tmpID,
                                    aCDKException);
                            tmpLogFileHandler.flush();
                            tmpSMILEScode = "[SMILES code could not be generated]";
                            tmpMinorExceptionsCounter++;
                        }
                        tmpOutput = tmpOutput.concat(SugarRemovalUtilityCmdApplication.OUTPUT_FILE_SEPARATOR + tmpSMILEScode);
                    }
                }
                tmpOutputFilePrintWriter.println(tmpOutput);
                tmpOutputFilePrintWriter.flush();
                tmpMoleculeCounter++;
            } catch (Exception anException) {
                SugarRemovalUtilityCmdApplication.LOGGER.log(Level.SEVERE,
                        anException.toString() + " Molecule number: " + tmpMoleculeCounter,
                        anException);
                tmpLogFileHandler.flush();
                tmpOutputFilePrintWriter.println(tmpOutput);
                tmpOutputFilePrintWriter.flush();
                tmpFatalExceptionCounter++;
                continue;
            }
        }
        System.out.println("Finished analysis of molecules.");
        System.out.println((tmpMoleculeCounter - 1) + " molecules were analysed.");
        System.out.println(tmpFatalExceptionCounter + " fatal exceptions occurred (molecule could not be read, was " +
                "disconnected or could not be cloned).");
        System.out.println(tmpMinorExceptionsCounter + " minor exceptions occurred (SMILES codes could not be generated etc.).");
        System.out.println(tmpSugarContainingMoleculesCounter + " molecules contained sugar moieties.");
        System.out.println(tmpContainsCircularSugarsCounter + " molecules contained circular sugar moieties.");
        System.out.println(tmpContainsLinearSugarsCounter + " molecules contained linear sugar moieties.");
        System.out.println(tmpBasicallyASugarCounter + " molecules were basically sugars because they were completely removed.");
        tmpFileInputStream.close();
        tmpOutputFilePrintWriter.flush();
        tmpOutputFilePrintWriter.close();
        tmpOutputFileWriter.close();
        tmpFileReader.close();
        tmpReader.close();
        tmpLogFileHandler.flush();
        tmpLogFileHandler.close();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Checks whether the given integer number is mapped to a type of moieties to remove, i.e. whether it is an accepted
     * input value for this parameter.
     *
     * @param aTypeOfMoietiesToRemove an integer that is supposed to be mapped to a type of moiety to remove
     * @return true if the given value is valid for the respective parameter
     */
    public static boolean isLegalTypeOfMoietiesToRemove(int aTypeOfMoietiesToRemove) {
        if (aTypeOfMoietiesToRemove == 1 || aTypeOfMoietiesToRemove == 2 || aTypeOfMoietiesToRemove == 3) {
            return true;
        } else {
            return false;
        }
    }
    //</editor-fold>
}
