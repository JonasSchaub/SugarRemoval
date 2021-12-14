/*
 * MIT License
 *
 * Copyright (c) 2021 Jonas Schaub, Achim Zielesny, Christoph Steinbeck, Maria Sorokina
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

/**
 * Main entry point for the Sugar Removal Utility command-line application. The main() method instantiates the
 * SugarRemovalUtilityCmdApplication class, calls its execute() function and measures the time it takes for execution.
 *
 * @author Jonas Schaub
 * @version 1.3.0.0
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
