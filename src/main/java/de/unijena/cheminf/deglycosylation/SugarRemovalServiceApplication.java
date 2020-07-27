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
import org.openscience.cdk.io.*;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.IChemFormatMatcher;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.*;
import java.util.List;
import java.util.Objects;
import java.util.logging.*;

/**
 * TODO
 */
public class SugarRemovalServiceApplication {

    /**
     * TODO
     */
    private SugarRemovalUtility sugarRemovalUtil;

    /**
     * TODO
     */
    private int typeOfMoietiesToRemove;

    /**
     * TODO
     */
    private File inputFile;

    /**
     *
     */
    private static final String OUTPUT_FILE_SEPARATOR = ";";

    /**
     *
     */
    private static final Logger LOGGER = Logger.getLogger(SugarRemovalServiceApplication.class.getName());

    /**
     * TODO
     */
    public SugarRemovalServiceApplication(File aFile, int aTypeOfMoietiesToRemove) throws IllegalArgumentException {
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
     * TODO
     */
    public SugarRemovalServiceApplication(File aFile,
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
        if (!SugarRemovalServiceApplication.isLegalTypeOfMoietiesToRemove(aTypeOfMoietiesToRemove)) {
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
        this.sugarRemovalUtil.setLinearSugarCandidateMinSizeSetting(aLinearSugarCandidateMinSizeSetting);
        this.sugarRemovalUtil.setLinearSugarCandidateMaxSizeSetting(aLinearSugarCandidateMaxSizeSetting);
        this.sugarRemovalUtil.setDetectLinearAcidicSugarsSetting(aDetectLinearAcidicSugarsSetting);
        this.sugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(aDetectSpiroRingsAsCircularSugarsSetting);
        this.sugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
    }

    /**
     * TODO
     */
    public static boolean isLegalTypeOfMoietiesToRemove(int aTypeOfMoietiesToRemove) {
        if (aTypeOfMoietiesToRemove == 1 || aTypeOfMoietiesToRemove == 2 || aTypeOfMoietiesToRemove == 3) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * TODO doc
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
        String tmpOutputFileHeader = "Nr" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "ID" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "originalMoleculeSMILES" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "deglycosylatedMoleculeSMILES" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "hadOrHasSugars" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "SugarMoietySMILES";
        tmpOutputFilePrintWriter.println(tmpOutputFileHeader);
        tmpOutputFilePrintWriter.flush();
        System.out.println("Results will be written to: " + tmpOutputFile.getAbsolutePath());
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputRootDirectoryPath + "Log.txt");
        } catch (IOException anIOException) {
            SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
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
        System.out.println("Detect circular sugars only with O-glycosidic bond: " + this.sugarRemovalUtil.areOnlyCircularSugarsWithOGlycosidicBondDetected());
        System.out.println("Remove only terminal sugar moieties: " + this.sugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println("Structure to keep mode setting: " + this.sugarRemovalUtil.getPreservationModeSetting().name()
                + " (" + this.sugarRemovalUtil.getPreservationModeSetting().ordinal() + ")");
        System.out.println("Structure to keep mode threshold: " + this.sugarRemovalUtil.getPreservationModeThresholdSetting());
        System.out.println("Detect circular sugars only with enough exocyclic oxygen atoms: " + this.sugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected());
        System.out.println("Exocyclic oxygen atoms to atoms in ring ratio threshold: " + this.sugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        System.out.println("Detect linear sugars in rings: " + this.sugarRemovalUtil.areLinearSugarsInRingsDetected());
        System.out.println("Linear sugar candidate minimum size: " + this.sugarRemovalUtil.getLinearSugarCandidateMinSizeSetting() + " carbon atoms");
        System.out.println("Linear sugar candidate maximum size: " + this.sugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting() + " carbon atoms");
        System.out.println("Detect linear acidic sugars: " + this.sugarRemovalUtil.areLinearAcidicSugarsDetected());
        System.out.println("Detect spiro rings as circular sugars: " + this.sugarRemovalUtil.areSpiroRingsDetectedAsCircularSugars());
        System.out.println("Entering iteration of molecules...");
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
                tmpOutput = tmpMoleculeCounter + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR;
                try {
                    tmpMolecule = tmpReader.next();
                } catch (Exception anException) {
                    SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE,
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
                tmpOutput = tmpOutput.concat(tmpID + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR);
                try {
                    tmpOriginalMoleculeSMILES = tmpSmiGen.create(tmpMolecule);
                } catch (CDKException aCDKException) {
                    SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE,
                            aCDKException.toString() + " Molecule id: " + tmpID,
                            aCDKException);
                    tmpLogFileHandler.flush();
                    tmpOriginalMoleculeSMILES = "[SMILES code could not be generated]";
                    tmpMinorExceptionsCounter++;
                }
                tmpOutput = tmpOutput.concat(tmpOriginalMoleculeSMILES + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR);
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
                    SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE,
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
                        SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE,
                                aCDKException.toString() + " Molecule id: " + tmpID,
                                aCDKException);
                        tmpLogFileHandler.flush();
                        tmpDeglycosylatedMoleculeSMILES = "[SMILES code could not be generated]";
                        tmpMinorExceptionsCounter++;
                    }
                }
                tmpOutput = tmpOutput.concat(tmpDeglycosylatedMoleculeSMILES + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                        + tmpHadSugars);
                if (!tmpSugarMoieties.isEmpty()) {
                    for (IAtomContainer tmpMoiety : tmpSugarMoieties) {
                        String tmpSMILEScode = null;
                        try {
                            tmpSMILEScode = tmpSmiGen.create(tmpMoiety);
                        } catch (CDKException aCDKException) {
                            SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE,
                                    aCDKException.toString() + " Molecule id: " + tmpID,
                                    aCDKException);
                            tmpLogFileHandler.flush();
                            tmpSMILEScode = "[SMILES code could not be generated]";
                            tmpMinorExceptionsCounter++;
                        }
                        tmpOutput = tmpOutput.concat(SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR + tmpSMILEScode);
                    }
                }
                tmpOutputFilePrintWriter.println(tmpOutput);
                tmpOutputFilePrintWriter.flush();
                tmpMoleculeCounter++;
            } catch (Exception anException) {
                SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE,
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
        System.out.println(tmpFatalExceptionCounter + " fatal exceptions occurred (molecule could not be read, was disconnected or could not be cloned.");
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
}
