/*
 * MIT License
 *
 * Copyright (c) 2020 Maria Sorokina, Jonas Schaub, Christoph Steinbeck
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

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.*;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.*;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

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
        this(aFile,
                aTypeOfMoietiesToRemove,
                SugarRemovalUtility.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_DEFAULT,
                SugarRemovalUtility.REMOVE_ONLY_TERMINAL_SUGARS_DEFAULT,
                SugarRemovalUtility.STRUCTURE_TO_KEEP_MODE_DEFAULT.ordinal(),
                SugarRemovalUtility.STRUCTURE_TO_KEEP_MODE_DEFAULT.getDefaultThreshold(),
                SugarRemovalUtility.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_DEFAULT,
                SugarRemovalUtility.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_DEFAULT,
                SugarRemovalUtility.DETECT_LINEAR_SUGARS_IN_RINGS_DEFAULT,
                SugarRemovalUtility.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_DEFAULT,
                SugarRemovalUtility.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_DEFAULT,
                SugarRemovalUtility.DETECT_LINEAR_ACIDIC_SUGARS_DEFAULT,
                SugarRemovalUtility.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_DEFAULT);
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
        SugarRemovalUtility.StructureToKeepModeOption tmpStructureToKeepMode = null;
        for (SugarRemovalUtility.StructureToKeepModeOption tmpEnumConstant : SugarRemovalUtility.StructureToKeepModeOption.values()) {
            if (aStructureToKeepModeSetting == tmpEnumConstant.ordinal()) {
                tmpStructureToKeepMode = tmpEnumConstant;
            }
        }
        if (Objects.isNull(tmpStructureToKeepMode)) {
            throw new IllegalArgumentException("Given structure to keep mode setting does not correspond to an ordinal " +
                    "in the StructureToKeepModeOption enumeration.");
        }
        this.sugarRemovalUtil.setStructureToKeepModeSetting(tmpStructureToKeepMode);
        if (aStructureToKeepModeThresholdSetting < 0) {
            throw new IllegalArgumentException("The structure to keep mode setting threshold " +
                    "cannot be negative.");
        }
        //case 1: Structure to keep mode setting is 'all' (ordinal 0). Then, a 0 threshold should be passed.
        // case 2: Structure to keep mode is not 'all' (ordinal nonzero). Then, a nonzero threshold should be passed.
        if (tmpStructureToKeepMode == SugarRemovalUtility.StructureToKeepModeOption.ALL && aStructureToKeepModeThresholdSetting != 0) {
            throw new IllegalArgumentException("Structure to keep mode 'all' was selected. Therefore, passing a nonzero " +
                    "threshold makes no sense.");
        } else if (tmpStructureToKeepMode != SugarRemovalUtility.StructureToKeepModeOption.ALL && aStructureToKeepModeThresholdSetting == 0) {
            throw new IllegalArgumentException("Please select structure to keep mode 'all' if the associated " +
                    "threshold should be 0.");
        }
        this.sugarRemovalUtil.setStructureToKeepModeThresholdSetting(aStructureToKeepModeThresholdSetting);
        this.sugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(aDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting);
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
     * TODO
     */
    public void execute() throws IOException, SecurityException, IllegalArgumentException {
        FileInputStream tmpFileInputStream = new FileInputStream(this.inputFile);
        InputStreamReader tmpInputStreamReader = new InputStreamReader(tmpFileInputStream);
        BufferedReader tmpBuffReader = new BufferedReader(tmpInputStreamReader);
        FormatFactory tmpCDKFormatFactory = new FormatFactory();
        IChemFormat tmpFormat = tmpCDKFormatFactory.guessFormat(tmpBuffReader);
        String tmpFormatClassName = tmpFormat.getFormatName();
        boolean tmpReaderIsIterating = false;
        IChemObjectReader tmpReader = null;
        switch (tmpFormatClassName) {
            case "MDL Structure-data file":
                tmpReaderIsIterating = true;
                tmpReader = new IteratingSDFReader(new FileReader(this.inputFile), DefaultChemObjectBuilder.getInstance(), true);
                break;
            case "SMILES":
                tmpReaderIsIterating = true;
                tmpReader = new IteratingSMILESReader(new FileReader(this.inputFile), DefaultChemObjectBuilder.getInstance());
                break;
            case "MDL Molfile":
            case "MDL Molfile V2000":
                tmpReaderIsIterating = false;
                tmpReader = new MDLV2000Reader(new FileReader(this.inputFile), IChemObjectReader.Mode.RELAXED);
                break;
            case "MDL Mol/SDF V3000":
                tmpReaderIsIterating = false;
                tmpReader = new MDLV3000Reader(new FileReader(this.inputFile), IChemObjectReader.Mode.RELAXED);
                break;
            default:
                throw new IllegalArgumentException("Given file type " + tmpFormatClassName + " cannot be used");
        }
        String tmpOutputRootDirectoryPath = this.inputFile.getAbsoluteFile().getParent() + File.separator;
        String tmpOutputFilePath = tmpOutputRootDirectoryPath + "deglycosylation_results.txt";
        File tmpOutputFile = new File(tmpOutputFilePath);
        FileWriter tmpOutputFileWriter = new FileWriter(tmpOutputFile);
        PrintWriter tmpOutputFilePrintWriter = new PrintWriter(tmpOutputFileWriter, true);
        String tmpOutputFileHeader = "ID" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "originalMoleculeSMILES" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "deglycosylatedMoleculeSMILES" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "hadOrHasSugars" + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                + "SugarMoietySMILES";
        tmpOutputFilePrintWriter.println(tmpOutputFileHeader);
        tmpOutputFilePrintWriter.flush();
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        IAtomContainer tmpMolecule = new AtomContainer();
        String tmpID;
        String tmpOriginalMoleculeSMILES;
        boolean tmpHadSugars;
        String tmpDeglycosylatedMoleculeSMILES;
        List<IAtomContainer> tmpSugarMoieties;
        if (!tmpReaderIsIterating) {
            try {
                if (tmpReader instanceof MDLV2000Reader) {
                    tmpMolecule = ((MDLV2000Reader) tmpReader).read(new AtomContainer());
                } else if (tmpReader instanceof MDLV3000Reader) {
                    tmpMolecule = ((MDLV3000Reader) tmpReader).read(new AtomContainer());
                } else {
                    throw new IllegalArgumentException("Something went wrong at file type determination.");
                }
            } catch (CDKException aCDKException) {
                SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE, aCDKException.toString(), aCDKException);
                tmpOutputFilePrintWriter.println("[Molecule could not be read]");
                tmpOutputFilePrintWriter.flush();
                return;
            }
            if (!Objects.isNull(tmpMolecule.getProperty(CDKConstants.TITLE))) {
                tmpID = tmpMolecule.getProperty(CDKConstants.TITLE);
            } else {
                tmpID = "[No title available]";
            }
            try {
                tmpOriginalMoleculeSMILES = tmpSmiGen.create(tmpMolecule);
            } catch (CDKException aCDKException) {
                SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE, aCDKException.toString(), aCDKException);
                tmpOriginalMoleculeSMILES = "[SMILES code could not be generated]";
            }
            if (this.sugarRemovalUtil.areOnlyTerminalSugarsRemoved()) {
                boolean tmpIsConnected = ConnectivityChecker.isConnected(tmpMolecule);
                if (!tmpIsConnected) {
                    tmpDeglycosylatedMoleculeSMILES = "[Deglycosylation not possible because structure is unconnected]";
                    tmpOutputFilePrintWriter.println(tmpID + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                            + tmpOriginalMoleculeSMILES + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                            + tmpDeglycosylatedMoleculeSMILES);
                    tmpOutputFilePrintWriter.flush();
                    return;
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
                SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE, aCloneNotSupportedException.toString(), aCloneNotSupportedException);
                tmpDeglycosylatedMoleculeSMILES = "[Deglycosylation not possible because structure could not be cloned]";
                tmpOutputFilePrintWriter.println(tmpID + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                        + tmpOriginalMoleculeSMILES + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                        + tmpDeglycosylatedMoleculeSMILES);
                tmpOutputFilePrintWriter.flush();
                return;
            }
            tmpHadSugars = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY);
            IAtomContainer tmpDeglycosylatedCore = tmpSugarMoieties.remove(0);
            if (tmpDeglycosylatedCore.isEmpty()) {
                tmpDeglycosylatedMoleculeSMILES = "";
            } else {
                try {
                    tmpDeglycosylatedMoleculeSMILES = tmpSmiGen.create(tmpDeglycosylatedCore);
                } catch (CDKException aCDKException) {
                    SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE, aCDKException.toString(), aCDKException);
                    tmpDeglycosylatedMoleculeSMILES = "[SMILES code could not be generated]";
                }
            }
            String tmpOutput = tmpID + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                    + tmpOriginalMoleculeSMILES + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                    + tmpDeglycosylatedMoleculeSMILES + SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR
                    + tmpHadSugars;
            if (!tmpSugarMoieties.isEmpty()) {
                for (IAtomContainer tmpMoiety : tmpSugarMoieties) {
                    String tmpSMILEScode = null;
                    try {
                        tmpSMILEScode = tmpSmiGen.create(tmpMoiety);
                    } catch (CDKException aCDKException) {
                        SugarRemovalServiceApplication.LOGGER.log(Level.SEVERE, aCDKException.toString(), aCDKException);
                        tmpSMILEScode = "[SMILES code could not be generated]";
                    }
                    tmpOutput = tmpOutput.concat(SugarRemovalServiceApplication.OUTPUT_FILE_SEPARATOR + tmpSMILEScode);
                }
            }
            tmpOutputFilePrintWriter.println(tmpOutput);
            tmpOutputFilePrintWriter.flush();
        } else {
            //TODO
        }
        tmpFileInputStream.close();
        tmpOutputFilePrintWriter.flush();
        tmpOutputFilePrintWriter.close();
        tmpOutputFileWriter.close();
        tmpReader.close();
    }
}
