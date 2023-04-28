/*
 * MIT License
 *
 * Copyright (c) 2023 Jonas Schaub, Achim Zielesny, Christoph Steinbeck, Maria Sorokina
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

import java.io.IOException;
import java.util.Locale;
import java.util.Objects;

/**
 * Main entry point for the Sugar Removal Utility command-line application. The main() method instantiates the
 * SugarRemovalUtilityCmdApplication class, calls its execute() function and measures the time it takes for execution.
 *
 * @author Jonas Schaub
 * @version 1.3.2.0
 */
public class Main {
    /**
     * Instantiates the SugarRemovalUtilityCmdApplication class and calls its execute() method.
     * This method also tests the Java version (must be version 11.0.5 or higher) and measures the time it takes to
     * execute the application. The command-line arguments (args) are not parsed or tested here but passed on to the
     * SugarRemovalUtilityCmdApplication class constructor. Outputs are written to System.out and System.err if
     * anything goes wrong. Exceptions are not thrown. The system is exited if an error occurs or the application is
     * run successfully.
     *
     * @param args the command line arguments, see SugarRemovalUtilityCmdApplication class constructor for more info
     */
    public static void main(String[] args) {
        try {
            System.out.println();
            String tmpJavaVersion = System.getProperty("java.version");
            if (Main.compareVersions(tmpJavaVersion, "17.0.4") < 0) {
                System.err.println("The version of your Java installation has to be at least 17.0.4 for this application to run.");
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
    //
    /**
     * Compares two version strings of the form "x.y.z" of variable length position by position after splitting at each
     * "." and parsing to integers. The result is analogous to what Integer.compare(int a, int b) returns for the
     * first position that differs between the two version strings. If they are of different size and the shorter
     * one is a substring of the longer one starting at position 0, their lengths are compared.
     *
     * @param aVersionString1 one version string v1
     * @param aVersionString2 another version string v2
     * @return the value 0 if v1 == v2; a value less than 0 if v1 smaller than v2; and a value greater than 0 if v1 greater than v2
     * @throws IllegalArgumentException if one of the parameters is null, empty, or blank
     * @author Jonas Schaub
     */
    public static int compareVersions(String aVersionString1, String aVersionString2) throws IllegalArgumentException {
        if (Objects.isNull(aVersionString1) || aVersionString1.isEmpty() || aVersionString1.isBlank()
                || Objects.isNull(aVersionString2) || aVersionString2.isEmpty() || aVersionString2.isBlank()) {
            throw new IllegalArgumentException("One of the arguments is null, empty or blank.");
        }
        String[] tmpSeparateNumbersV1 = aVersionString1.split("\\.");
        String[] tmpSeparateNumbersV2 = aVersionString2.split("\\.");
        int tmpIterations = tmpSeparateNumbersV1.length < tmpSeparateNumbersV2.length ? tmpSeparateNumbersV1.length : tmpSeparateNumbersV2.length;
        for (int i = 0; i < tmpIterations; i++) {
            int tmpV1Int = Integer.parseInt(tmpSeparateNumbersV1[i]);
            int tmpV2Int = Integer.parseInt(tmpSeparateNumbersV2[i]);
            int tmpResult = Integer.compare(tmpV1Int, tmpV2Int);
            if (tmpResult != 0) {
                return tmpResult;
            }
        }
        return Integer.compare(tmpSeparateNumbersV1.length, tmpSeparateNumbersV2.length);
    }
}
