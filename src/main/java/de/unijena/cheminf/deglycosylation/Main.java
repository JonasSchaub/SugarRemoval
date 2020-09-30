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
 * - update doc
 * - update version
 */

import java.io.IOException;
import java.util.Locale;

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
            System.out.println();
            String tmpJavaVersion = System.getProperty("java.version");
            if (tmpJavaVersion.compareTo("11.0.5") < 0) {
                System.err.println("The version of your Java installation has to be at least 11.0.5 for this application to run.");
                System.exit(-1);
            }
            Locale.setDefault(Locale.US);
            System.out.println("Sugar Removal Service Application starting. Evaluating command line arguments...");
            SugarRemovalUtilityCmdApplication tmpSugarRemovalApp = null;
            try {
                tmpSugarRemovalApp = new SugarRemovalUtilityCmdApplication(args);
            } catch (IllegalArgumentException anIllegalArgumentException) {
                System.err.println("An error occurred while parsing the command-line arguments: " + anIllegalArgumentException.getMessage());
                System.err.println("Use the -h or --help option to print usage and help information. Use -v or --version " +
                        "to print the version of your SRU CMD App JAR.");
                System.err.println("Application will now exit.");
                System.exit(-1);
            }
            if (tmpSugarRemovalApp.wasHelpOrVersionQueried()) {
                System.exit(0);
            }
            System.out.println("All command line arguments valid. Executing application...");
            long tmpStartTime = System.currentTimeMillis();
            try {
                tmpSugarRemovalApp.execute();
            } catch (IOException | SecurityException | IllegalArgumentException | NullPointerException anException) {
                System.err.println("Problem at execution of application: " + anException.getMessage());
                System.exit(-1);
            }
            long tmpEndTime = System.currentTimeMillis();
            double tmpDurationInMinutes = ((double)tmpEndTime - (double)tmpStartTime) / 60000;
            System.out.println("Execution successful. Analysis took " + String.format("%.3f", tmpDurationInMinutes) + " minutes.");
            System.out.println("Application will now exit.");
            System.exit(0);
        } catch (Exception anException) {
            System.err.println("An unexpected error occurred: " + anException.getMessage());
            anException.printStackTrace(System.err);
            System.exit(-1);
        }
    }
}
