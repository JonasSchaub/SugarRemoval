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

import java.io.File;
import java.util.Objects;

/**
 * TODO add doc
 */
public class Main {
    /**
     * TODO
     * See Readers in COCONUT
     * Output: Meditate about it (the detailed structure of the output file) CSV?
     * Should include the SMILES of the original structure and the removed moieties
     *
     * arg 0: Path to SMILES file, SDF, or molfile [String]
     * arg 1: param to indicate whether circular, linear, or both types of sugars should be removed [int?]
     * the other params should be optional (but if one is present, all should be present?)
     *
     * arg 2: detect circular sugars only with O glycosidic bond [boolean]
     * arg 3: remove only terminal sugars [boolean]
     *      -> this requires checking all the molecules for being connected before processing them!
     *      -> must be done by the app, not here
     * arg 4: structure to keep mode [int?]
     * arg 5: structure to keep mode threshold [int]
     * arg 6: detect circular sugars only with enough exocyclic oxygen atoms [boolean]
     * arg 7: ratio attached oxygen atoms to atoms in ring threshold [double]
     * arg 8: detect linear sugars in rings [boolean]
     * arg 9: linear sugar candidate min size [int]
     * arg 10: linear sugar candidate max size [int]
     * arg 11: detect linear acidic sugars [boolean]
     * arg 12: detect spiro rings as circular sugars [boolean]
     */
    public static void main(String[] args) {
        try {
            if (args.length != 2 && args.length != 13) {
                System.err.println("Number of command line arguments must be either 2 or 13.");
                System.exit(-1);
            }
            SugarRemovalServiceApplication tmpSugarRemovalApp = null;
            String tmpPath = args[0];
            if (Objects.isNull(tmpPath) || tmpPath.isBlank()) {
                System.err.println("Given path to SMILES, SD, or MOL file (argument at position 0) is empty or blank.");
            }
            File tmpFile = new File(tmpPath);
            if (!tmpFile.exists() || !tmpFile.canRead() || !tmpFile.isFile()) {
                System.err.println("Given path to SMILES, SD, or MOL file (argument at position 0) does not exist or " +
                        "the file cannot be read.");
                System.exit(-1);
            }
            int tmpTypeOfMoietiesToRemove = -1;
            try {
                tmpTypeOfMoietiesToRemove = Integer.parseInt(args[1]);
            } catch (NumberFormatException aNumberFormatException) {
                System.err.println("Integer indicating which type of moieties is to remove (argument at position 1)" +
                        "cannot be parsed.");
                System.exit(-1);
            }
            if (!SugarRemovalServiceApplication.isLegalTypeOfMoietiesToRemove(tmpTypeOfMoietiesToRemove)) {
                System.err.println("Argument at position 1 indicating which type of moieties is to remove must be " +
                        "either 1 (circular sugar moieties), 2 (linear sugar moieties), or 3 (both).");
                System.exit(-1);
            }
            if (args.length == 2) {
                tmpSugarRemovalApp = new SugarRemovalServiceApplication(tmpFile, tmpTypeOfMoietiesToRemove);
            } else if (args.length == 13) {
                //note: The boolean returned represents the value true if the string argument is not null and is equal,
                // ignoring case, to the string "true". Otherwise, a false value is returned, including for a null argument.
                // So, e.g. 'hugo' results in false and does not cause an exception
                boolean tmpDetectCircularSugarsOnlyWithOGlycosidicBondSetting = Boolean.parseBoolean(args[2]);
                boolean tmpRemoveOnlyTerminalSugarsSetting = Boolean.parseBoolean(args[3]);
                int tmpStructureToKeepModeSetting = -1;
                try {
                    tmpStructureToKeepModeSetting = Integer.parseInt(args[4]);
                } catch (NumberFormatException aNumberFormatException) {
                    System.err.println("Integer indicating the structure to keep mode setting (argument at position 4)" +
                            "cannot be parsed.");
                    System.exit(-1);
                }
                //the given int is taken as ordinal value of the indicated enum object; ordinals start at 0
                if (tmpStructureToKeepModeSetting >= SugarRemovalUtility.StructureToKeepModeOption.values().length ||
                        tmpStructureToKeepModeSetting < 0) {
                    System.err.println("Argument at position 4 indicating the structure to keep mode setting must be " +
                            "either 0 (keep all), 1 (judge by heavy atom count), or 2 (judge by molecular weight).");
                    System.exit(-1);
                }
                int tmpStructureToKeepModeThresholdSetting = -1;
                try {
                    tmpStructureToKeepModeThresholdSetting = Integer.parseInt(args[5]);
                } catch (NumberFormatException aNumberFormatException) {
                    System.err.println("Integer indicating the structure to keep mode setting threshold (argument at position 5)" +
                            "cannot be parsed.");
                    System.exit(-1);
                }
                if (tmpStructureToKeepModeThresholdSetting < 0) {
                    System.err.println("Integer indicating the structure to keep mode setting threshold (argument at position 5)" +
                            "cannot be negative.");
                    System.exit(-1);
                }
                //case 1: Structure to keep mode setting is 'all' (ordinal 0). Then, a 0 threshold should be passed.
                // case 2: Structure to keep mode is not 'all' (ordinal nonzero). Then, a nonzero threshold should be passed.
                if (tmpStructureToKeepModeSetting == 0 && tmpStructureToKeepModeThresholdSetting != 0) {
                    System.err.println("Structure to keep mode 'all' was selected (argument at position 4 is '0'). " +
                            "Therefore, passing a nonzero threshold (argument at position 5) makes no sense.");
                    System.exit(-1);
                } else if (tmpStructureToKeepModeSetting != 0 && tmpStructureToKeepModeThresholdSetting == 0) {
                    System.err.println("Please select structure to keep mode 'all' (argument value 0 at position 4) if " +
                            "the associated threshold should be 0 (argument at position 5.");
                    System.exit(-1);
                }
                boolean tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting = Boolean.parseBoolean(args[6]);
                double tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting = -1;
                try {
                    tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting = Double.parseDouble(args[7]);
                } catch (NumberFormatException | NullPointerException anException) {
                    System.err.println("Number indicating the exocyclic oxygen atoms to atoms in ring ratio threshold " +
                            "(argument at position 7) cannot be parsed.");
                    System.exit(-1);
                }
                boolean tmpIsFinite = Double.isFinite(tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting);
                boolean tmpIsNegative = (tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting < 0);
                if (!tmpIsFinite || tmpIsNegative) {
                    System.err.println("Number indicating the exocyclic oxygen atoms to atoms in ring ratio threshold" +
                            "(argument at position 7) is NaN, infinite or negative.");
                    System.exit(-1);
                }
                //case 1: If the number of exocyclic oxygen atoms is neglected, a 0 threshold should be passed
                // case 2: If it is detected, a nonzero threshold should be passed
                if (!tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
                        && tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting != 0) {
                    System.err.println("The number of exocyclic oxygen atoms of circular sugars is neglected at detection " +
                            "(argument at position 6 is 'false'). " +
                            "Therefore, passing a nonzero threshold (argument at position 7) makes no sense.");
                    System.exit(-1);
                } else if (tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
                        && tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting == 0) {
                    System.err.println("Please select to neglect the number of exocyclic oxygen atoms (argument value " +
                            "'false' at position 6) if the associated threshold should be 0 (argument at position 7).");
                    System.exit(-1);
                }
                boolean tmpDetectLinearSugarsInRingsSetting = Boolean.parseBoolean(args[8]);
                int tmpLinearSugarCandidateMinSizeSetting = -1;
                int tmpLinearSugarCandidateMaxSizeSetting = -1;
                try {
                    tmpLinearSugarCandidateMinSizeSetting = Integer.parseInt(args[9]);
                    tmpLinearSugarCandidateMaxSizeSetting = Integer.parseInt(args[10]);
                } catch (NumberFormatException aNumberFormatException) {
                    System.err.println("The linear sugar candidate minimum (argument at position 9) or maximum size " +
                            "(argument at position 10) cannot be parsed.");
                    System.exit(-1);
                }
                if (tmpLinearSugarCandidateMinSizeSetting < 1 || tmpLinearSugarCandidateMaxSizeSetting < 1) {
                    System.err.println("The linear sugar candidate minimum and maximum sizes must both be equal or " +
                            "higher than 1 (arguments at positions 9 and 10).");
                    System.exit(-1);
                }
                boolean tmpDetectLinearAcidicSugarsSetting = Boolean.parseBoolean(args[11]);
                boolean tmpDetectSpiroRingsAsCircularSugarsSetting = Boolean.parseBoolean(args[12]);
                tmpSugarRemovalApp = new SugarRemovalServiceApplication(tmpFile,
                        tmpTypeOfMoietiesToRemove,
                        tmpDetectCircularSugarsOnlyWithOGlycosidicBondSetting,
                        tmpRemoveOnlyTerminalSugarsSetting,
                        tmpStructureToKeepModeSetting,
                        tmpStructureToKeepModeThresholdSetting,
                        tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting,
                        tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting,
                        tmpDetectLinearSugarsInRingsSetting,
                        tmpLinearSugarCandidateMinSizeSetting,
                        tmpLinearSugarCandidateMaxSizeSetting,
                        tmpDetectLinearAcidicSugarsSetting,
                        tmpDetectSpiroRingsAsCircularSugarsSetting);
            } else {
                //redundant...
                System.err.println("Number of command line arguments must be either 2 or 13.");
                System.exit(-1);
            }
            tmpSugarRemovalApp.execute();
            System.exit(0);
        } catch (Exception anException) {
            System.err.println("An unexpected error occurred: ");
            anException.printStackTrace(System.err);
            System.exit(-1);
        }
    }
}
