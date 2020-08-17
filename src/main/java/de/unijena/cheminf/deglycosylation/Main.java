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

/**
 * TODO
 * - improve command line argument management
 */

import java.io.File;
import java.io.IOException;
import java.util.Locale;
import java.util.Objects;

/**
 * Main entry point for the SugarRemovalUtility command line application. The main() method parses the command line
 * arguments, instantiates, and starts the application.
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class Main {
    /**
     * Parses the command line arguments, instantiates the SugarRemovalUtilityCmdApplication, and executes it.
     * The method also tests the Java version (must be version 11.0.5 or higher) and measures the time it takes to
     * execute the application. Outputs are written to System.out and System.err if anything goes wrong. Exceptions
     * are not thrown. The system is exited if an error occurs or the application is run successfully.
     * <br>The command line arguments 'args' must be constructed as follows:
     * <p>
     * * args[0]: Path to the input file, either absolute or relative to the current directory. Example:
     * "D:\Project_Sugar_Removal\SugarRemovalUtility CMD App\smiles_test_file.txt" or
     * "smiles_test_file.txt" if the console is already in the "SugarRemovalUtility CMD App"
     * directory. The backslahes '\' are used in a Microsoft Windows operating system; it should be slash
     * '/' in a Unix shell. Double quotes " are not mandatory but recommended to allow for spaces in the path. The
     * path must not be empty and the given file must exist and be accessible and readable. The file type extension
     * is not important for the determination of the file type but it must be specified in the path. Accepted input
     * formats: MDL Molfile, MDL Structure data file (SDF) and SMILES files (of format: [SMILES string][space][name]
     * in each line, see example file).
     * <br>* args[1]: A number of ["1","2","3"] indicating which type of sugar moieties should be removed, "1" for
     * circular sugar moieties, "2" for linear sugar moieties, or "3" for circular AND linear sugar moieties.
     * <p>
     * These two arguments must ALWAYS be given. The remaining arguments are optional. But either
     * none of them or all of them must be given. If only these two arguments are given, the other
     * settings will be in their default value.
     * <p>
     * * args[2]: Either "true" or "false", indicating whether circular sugars should be detected (and
     * removed) only if they have an O-glycosidic bond to another moiety or the core of the molecule. Any other value
     * of this argument will be interpreted as "false". Default: "false".
     * <br>* args[3]: Either "true" or "false", indicating whether only terminal sugar moieties should be removed.
     * Any other value of this argument will be interpreted as "false". Default: "true". Important note: If this
     * setting is set to "true", the input molecules must all consist of one connected structure, respectively. If they
     * already contain multiple, disconnected structures (e.g. counter-ions), the respective molecules are ignored.
     * <br>* args[4]: A number of ["0","1","2"] indicating which preservation mode to use. This specifies under what
     * circumstances to discard structures that get disconnected from the central core in the sugar removal process,
     * "0" to preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected,
     * depending on the other settings), "1" to remove disconnected structures that do not have enough heavy atoms, or
     * "2" to remove disconnected structures that do not have a sufficient molecular weight. Default: "1" (judge
     * disconnected structures by their heavy atom count).
     * <br>* args[5]: An integer number giving the threshold of the preservation mode, i.e. how many heavy atoms a
     * disconnected structure needs to have at least to be not removed or how heavy (in terms of its molecular weight)
     * it needs to be. Default: "5" (heavy atoms). The integer number must be positive. If the previous argument at
     * position 4 was passed the value "0" (preserve all structures), this argument must also be passed a zero value.
     * In the opposite case, this argument must be passed a non-zero value if the previous argument at position 4 was
     * given the value 1 or 2.
     * <br>* args[6]: Either "true" or "false", indicating whether circular sugars should be detected (and
     * removed) only if they have a sufficient number of attached exocyclic oxygen atoms. Any other value of this
     * argument will be interpreted as "false". Default: "true".
     * <br>* args[7]: A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring
     * ratio a circular sugar needs to have to be detected as such.
     * Default: "0.5" (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms).
     * If the previous argument at position 6 was passed the value "false" (detect circular sugars
     * neglecting their number of attached exocyclic oxygen atoms), this argument must be passed a zero
     * value. In the opposite case, this argument must be passed a non-zero value if the previous argument at
     * position 6 was given the value "true". The number must be positive.
     * <br>* args[8]: Either "true" or "false", indicating whether linear sugars in rings should be detected (and
     * removed). Any other value of this argument will be interpreted as "false". Default: "false".
     * <br>* args[9]: An integer number indicating the minimum number of carbon atoms a linear sugar needs to have
     * to be detected as such. Default: "4". The integer number must be positive and higher than or equal to 1 and also
     * smaller than the linear sugar candidate maximum size (argument at position 10).
     * <br>* args[10]: An integer number indicating the maximum number of carbon atoms a linear sugar needs to
     * have to be detected as such. Default: "7". The integer number must be positive and higher than or equal to 1 and
     * also higher than the linear sugar candidate minimum size (argument at previous position 9).
     * <br>* args[11]: Either "true" or "false", indicating whether linear acidic sugars should be included in the
     * set of linear sugar patterns for the initial detection. Any other value of this argument will be interpreted as
     * "false". Default: "false".
     * <br>* args[12]: Either "true" or "false", indicating whether spiro rings (rings that share one atom with
     * another cycle) should be included in the circular sugar detection. Any other value of this argument will be
     * interpreted as "false". Default: "false".
     * <p>
     * Example (all settings in default): String[] args = new String[] {"smiles_test_file.txt", "3", "false", "true", "1", "5", "true", "0.5",
     * "false", "4", "7", "false", "false"}
     *
     * @param args the command line arguments, see above
     */
    public static void main(String[] args) {
        try {
            String tmpJavaVersion = System.getProperty("java.version");
            if (tmpJavaVersion.compareTo("11.0.5") < 0) {
                System.err.println("The version of your Java installation has to be at least 11.0.5 for this application to run.");
                System.exit(-1);
            }
            System.out.println("Sugar Removal Service Application starting. Evaluating command line arguments...");
            if (args.length != 2 && args.length != 13) {
                System.err.println("Number of command line arguments must be either 2 or 13.");
                System.exit(-1);
            }
            Locale.setDefault(Locale.US);
            SugarRemovalUtilityCmdApplication tmpSugarRemovalApp = null;
            String tmpPath = args[0].trim();
            if (Objects.isNull(tmpPath) || tmpPath.isBlank()) {
                System.err.println("Given path to SMILES, SD, or MOL file (argument at position 0) is empty or blank.");
                System.exit(-1);
            }
            File tmpFile = new File(tmpPath);
            if (!tmpFile.exists() || !tmpFile.canRead() || !tmpFile.isFile()) {
                System.err.println("Given path to SMILES, SD, or MOL file (argument at position 0) does not exist or " +
                        "the file cannot be read.");
                System.exit(-1);
            }
            int tmpTypeOfMoietiesToRemove = -1;
            try {
                tmpTypeOfMoietiesToRemove = Integer.parseInt(args[1].trim());
            } catch (NumberFormatException aNumberFormatException) {
                System.err.println("Integer indicating which type of moieties is to remove (argument at position 1)" +
                        "cannot be parsed.");
                System.exit(-1);
            }
            if (!SugarRemovalUtilityCmdApplication.isLegalTypeOfMoietiesToRemove(tmpTypeOfMoietiesToRemove)) {
                System.err.println("Argument at position 1 indicating which type of moieties is to remove must be " +
                        "either 1 (circular sugar moieties), 2 (linear sugar moieties), or 3 (both).");
                System.exit(-1);
            }
            if (args.length == 2) {
                tmpSugarRemovalApp = new SugarRemovalUtilityCmdApplication(tmpFile, tmpTypeOfMoietiesToRemove);
            } else if (args.length == 13) {
                //note: The boolean returned represents the value true if the string argument is not null and is equal,
                // ignoring case, to the string "true". Otherwise, a false value is returned, including for a null argument.
                // So, e.g. 'hugo' results in false and does not cause an exception
                boolean tmpDetectCircularSugarsOnlyWithOGlycosidicBondSetting = Boolean.parseBoolean(args[2].trim());
                boolean tmpRemoveOnlyTerminalSugarsSetting = Boolean.parseBoolean(args[3].trim());
                int tmpStructureToKeepModeSetting = -1;
                try {
                    tmpStructureToKeepModeSetting = Integer.parseInt(args[4].trim());
                } catch (NumberFormatException aNumberFormatException) {
                    System.err.println("Integer indicating the structure to keep mode setting (argument at position 4)" +
                            "cannot be parsed.");
                    System.exit(-1);
                }
                //the given int is taken as ordinal value of the indicated enum object; ordinals start at 0
                if (tmpStructureToKeepModeSetting >= SugarRemovalUtility.PreservationModeOption.values().length ||
                        tmpStructureToKeepModeSetting < 0) {
                    System.err.println("Argument at position 4 indicating the structure to keep mode setting must be " +
                            "either 0 (keep all), 1 (judge by heavy atom count), or 2 (judge by molecular weight).");
                    System.exit(-1);
                }
                int tmpStructureToKeepModeThresholdSetting = -1;
                try {
                    tmpStructureToKeepModeThresholdSetting = Integer.parseInt(args[5].trim());
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
                boolean tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting = Boolean.parseBoolean(args[6].trim());
                double tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting = -1;
                try {
                    tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting = Double.parseDouble(args[7].trim());
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
                boolean tmpDetectLinearSugarsInRingsSetting = Boolean.parseBoolean(args[8].trim());
                int tmpLinearSugarCandidateMinSizeSetting = -1;
                int tmpLinearSugarCandidateMaxSizeSetting = -1;
                try {
                    tmpLinearSugarCandidateMinSizeSetting = Integer.parseInt(args[9].trim());
                    tmpLinearSugarCandidateMaxSizeSetting = Integer.parseInt(args[10].trim());
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
                if (tmpLinearSugarCandidateMinSizeSetting > tmpLinearSugarCandidateMaxSizeSetting) {
                    System.err.println("The linear sugar candidate minimum size (argument at position 9) must be smaller " +
                            "than the maximum size (argument at position 10).");
                    System.exit(-1);
                }
                boolean tmpDetectLinearAcidicSugarsSetting = Boolean.parseBoolean(args[11].trim());
                boolean tmpDetectSpiroRingsAsCircularSugarsSetting = Boolean.parseBoolean(args[12].trim());
                //throws IllegalArgumentException but that should not happen since all arguments have been thoroughly tested.
                tmpSugarRemovalApp = new SugarRemovalUtilityCmdApplication(tmpFile,
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
                System.out.println("All command line arguments valid. Executing application...");
            } else {
                //redundant...
                System.err.println("Number of command line arguments must be either 2 or 13.");
                System.exit(-1);
            }
            long tmpStartTime = System.currentTimeMillis();
            try {
                tmpSugarRemovalApp.execute();
            } catch (IOException | SecurityException | IllegalArgumentException anException) {
                System.err.println("Problem at execution of application: " + anException.toString());
                System.exit(-1);
            }
            long tmpEndTime = System.currentTimeMillis();
            double tmpDurationInMinutes = ((double)tmpEndTime - (double)tmpStartTime) / 60000;
            System.out.println("Execution successful. Analysis took " + String.format("%.3f", tmpDurationInMinutes) + " minutes.");
            System.out.println("Application will now exit.");
            System.exit(0);
        } catch (Exception anException) {
            System.err.println("An unexpected error occurred: ");
            anException.printStackTrace(System.err);
            System.exit(-1);
        }
    }
}
