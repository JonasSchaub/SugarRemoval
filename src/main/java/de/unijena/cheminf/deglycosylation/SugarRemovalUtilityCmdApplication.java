/*
 * MIT License
 *
 * Copyright (c) 2022 Jonas Schaub, Achim Zielesny, Christoph Steinbeck, Maria Sorokina
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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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
 * Controller of the Sugar Removal Utility command-line application. It can be used to remove sugar moieties from
 * molecules in a given data set, according to <a href="https://doi.org/10.1186/s13321-020-00467-y"
 * >"Schaub, J., Zielesny, A., Steinbeck, C., Sorokina, M. Too sweet: cheminformatics for deglycosylation in natural
 * products. J Cheminform 12, 67 (2020). https://doi.org/10.1186/s13321-020-00467-y"</a>. This class
 * basically instantiates the SugarRemovalUtility class with the settings specified in the command line arguments and
 * uses it to iterate over all molecules found in the given file and remove their sugar moieties. Also, a CSV file detailing
 * the deglycosylated cores and removed sugar moieties for each molecule is created as output.
 *
 * @author Jonas Schaub
 * @version 1.3.0.0
 */
public class SugarRemovalUtilityCmdApplication {
    //<editor-fold desc="Public static final constants">
    /**
     * Version string of this class to print out if -v --version is queried from the command-line.
     */
    public static final String VERSION = "1.3.0.0";
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

    //<editor-fold desc="Command-line options">
    /**
     * Collection of options for first stage of CMD argument parsing. If the usage instructions (-h --help) or version
     * of this class (-v --version) where queried, they have to be parsed without regarding the required options in
     * the second stage.
     */
    private static final Options HELP_AND_VERSION_OPTIONS = new Options();

    /**
     * Collection of options for second/main stage of CMD argument parsing. I.e. the 'true' CMD arguments of the
     * application (excluding help or version query).
     */
    private static final Options CMD_OPTIONS = new Options();

    /**
     * Command-line option to query usage information.
     */
    private static final Option HELP_OPTION = Option.builder("h")
            .longOpt("help")
            .hasArg(false)
            .required(false)
            .optionalArg(false)
            .desc("Print options and descriptions")
            .build();

    /**
     * Command-line option to query the version string of the application.
     */
    private static final Option VERSION_OPTION = Option.builder("v")
            .longOpt("version")
            .hasArg(false)
            .required(false)
            .optionalArg(false)
            .desc("Print version of the application")
            .build();

    /**
     * Required command-line option to pass the input file path.
     */
    private static final Option INPUT_FILE_PATH_OPTION = Option.builder("i")
            .argName("filePath")
            .longOpt("inputFilePath")
            .desc("Path to the input file; either absolute or relative to the current directory; required option")
            .hasArg(true)
            .numberOfArgs(1)
            .required(true)
            .optionalArg(false)
            .build();

    /**
     * Required command-line option to specify the type of moieties to remove.
     */
    private static final Option TYPE_OF_MOIETIES_TO_REMOVE_OPTION = Option.builder("t")
            .argName("integer")
            .longOpt("typeOfMoietiesToRemove")
            .desc("A number of [1 (circular), 2 (linear), 3 (both)] indicating which type of sugar moieties should be " +
                    "removed; required option")
            .hasArg(true)
            .numberOfArgs(1)
            .required(true)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to detect only circular sugars with an O-glycosidic bond.
     */
    private static final Option DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION = Option.builder("glyBond")
            .argName("boolean")
            .longOpt("detectCircSugOnlyWithGlyBond")
            .desc("Either true or false, indicating whether circular sugars should be detected (and " +
                    "removed) only if they have an O-glycosidic bond to another moiety or the core of the molecule; " +
                    "default: false")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to remove only terminal sugar moieties.
     */
    private static final Option REMOVE_ONLY_TERMINAL_SUGARS_OPTION = Option.builder("remTerm")
            .argName("boolean")
            .longOpt("removeOnlyTerminalSugars")
            .desc("Either true or false, indicating whether only terminal sugar moieties should be removed; " +
                    "default: true")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to set the preservation mode setting.
     */
    private static final Option PRESERVATION_MODE_OPTION = Option.builder("presMode")
            .argName("integer")
            .longOpt("preservationModeOption")
            .desc("A number of [0 (preserve all), 1 (judge by heavy atom count), 2 (judge by molecular weight)] indicating " +
                    "which preservation mode to use. This specifies under what circumstances to preserve or discard " +
                    "structures that get disconnected from the central core in the sugar removal process; " +
                    "default: 1 (judge by heavy atom count)")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to set the preservation mode threshold.
     */
    private static final Option PRESERVATION_MODE_THRESHOLD_OPTION = Option.builder("presThres")
            .argName("integer")
            .longOpt("preservationModeThreshold")
            .desc("An integer number giving the threshold of the preservation mode, i.e. how many heavy atoms a " +
                    "disconnected structure needs to have at least to be not " +
                    "removed or how heavy (in terms of its molecular weight) it needs to be; " +
                    "default: 5 (heavy atoms)")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to detect only circular sugars with a sufficient number of
     * exocyclic oxygen atoms attached to the central ring.
     */
    private static final Option DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION = Option.builder("oxyAtoms")
            .argName("boolean")
            .longOpt("detectCircSugOnlyWithEnoughExocycOxyAtoms")
            .desc("Either true or false, indicating whether circular sugars should be detected (and " +
                    "removed) only if they have a sufficient number of attached exocyclic oxygen atoms; " +
                    "default: true")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify the minimum ratio of attached exocyclic oxygen atoms to the number of
     * atoms in the possible sugar ring for circular sugars.
     */
    private static final Option EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION = Option.builder("oxyAtomsThres")
            .argName("number")
            .longOpt("exocycOxyAtomsToAtomsInRingRatioThreshold")
            .desc("A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring " +
                    "ratio a circular sugar needs to have to be detected as such; " +
                    "default: 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms)")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to detect linear sugars in rings.
     */
    private static final Option DETECT_LINEAR_SUGARS_IN_RINGS_OPTION = Option.builder("linSugInRings")
            .argName("boolean")
            .longOpt("detectLinSugInRings")
            .desc("Either true or false, indicating whether linear sugars in rings should be detected (and " +
                    "removed); " +
                    "default: false")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify the minimum size (number of carbon atoms) of linear sugars.
     */
    private static final Option LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION = Option.builder("linSugMinSize")
            .argName("integer")
            .longOpt("linSugCandidateMinimumSize")
            .desc("An integer number indicating the minimum number of carbon atoms a linear sugar needs to have " +
                    "to be detected as such; " +
                    "default: 4")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify the maximum size (number of carbon atoms) of linear sugars.
     */
    private static final Option LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION = Option.builder("linSugMaxSize")
            .argName("integer")
            .longOpt("linSugCandidateMaximumSize")
            .desc("An integer number indicating the Maximum number of carbon atoms a linear sugar needs to have " +
                    "to be detected as such; " +
                    "default: 7")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to add the linear acidic sugars to the linear sugar pattern set
     * for initial detection.
     */
    private static final Option DETECT_LINEAR_ACIDIC_SUGARS_OPTION = Option.builder("linAcSug")
            .argName("boolean")
            .longOpt("detectLinAcidicSug")
            .desc("Either true or false, indicating whether linear acidic sugars should be included in the " +
                    "set of linear sugar patterns for the initial detection; " +
                    "default: false")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to detect spiro rings as circular sugars.
     */
    private static final Option DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION = Option.builder("circSugSpiro")
            .argName("boolean")
            .longOpt("detectSpiroRingsAsCircSug")
            .desc("Either true or false, indicating whether spiro rings (rings that share one atom with " +
                    "another cycle) should be included in the circular sugar detection; " +
                    "default: false")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();

    /**
     * Optional command-line option to specify whether to detect circular sugar-like moieties having keto groups.
     */
    private static final Option DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION = Option.builder("circSugKetoGroups")
            .argName("boolean")
            .longOpt("detectCircularSugarsWithKetoGroups")
            .desc("Either true or false, indicating whether circular sugar-like moieties with keto groups should be " +
                    "detected; default: false")
            .hasArg(true)
            .numberOfArgs(1)
            .required(false)
            .optionalArg(false)
            .build();
    //</editor-fold>
    //</editor-fold>
    //
    //<editor-fold desc="Static initializer">
    /**
     * Static initializer that adds all the command-line options to the respective Options containers.
     */
    static {
        SugarRemovalUtilityCmdApplication.HELP_AND_VERSION_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.HELP_OPTION);
        SugarRemovalUtilityCmdApplication.HELP_AND_VERSION_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.VERSION_OPTION);

        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.REMOVE_ONLY_TERMINAL_SUGARS_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.DETECT_LINEAR_SUGARS_IN_RINGS_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.DETECT_LINEAR_ACIDIC_SUGARS_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION);
        SugarRemovalUtilityCmdApplication.CMD_OPTIONS.addOption(
                SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    /**
     * Boolean property specifying whether the usage instructions or version string were queried in the command-line
     * arguments. If this is the case, the application cannot be executed but should be exited.
     */
    private boolean wasHelpOrVersionQueried;

    /**
     * SugarRemovalUtility instance used to detect and remove the sugar moieties.
     */
    private SugarRemovalUtility sugarRemovalUtil;

    /**
     * Setting that specifies whether to remove circular, linear, or both types of moieties from the molecules.
     */
    private int typeOfMoietiesToRemove;

    /**
     * Input file containing molecule set, either MDL Molfile, MDL Structure data file (SDF), or SMILES file.
     */
    private File inputFile;
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    /**
     * Constructor that parses the command-line arguments and instantiates this class and the SugarRemovalUtility with
     * the given settings. If the -h --help or -v --version options are given, the remaining arguments are not parsed
     * and the application should be exited.
     * <br>The command line arguments 'args' must be constructed as follows:
     * <p>
     * * option -h --help: Print usage and help information regarding the required command-line arguments and options.
     * If this option is used, the constructor is exited afterwards.
     * <br>* option -v --version: Print version string of the Sugar Removal Utility Command-Line Application. If this
     * option is used, the constructor is exited afterwards.
     * <p>
     * * option -i --inputFilePath <filePath>: Path to the input file, either absolute or relative to the current
     * directory. Example: "D:\Project_Sugar_Removal\SugarRemovalUtility CMD App\smiles_test_file.txt" or
     * "smiles_test_file.txt" if the console is already in the "SugarRemovalUtility CMD App"
     * directory. The backslahes '\' are used in a Microsoft Windows operating system; it should be slash
     * '/' in a Unix shell. Double quotes " are not mandatory but recommended to allow for spaces in the path. The
     * path must not be empty and the given file must exist and be accessible and readable. The file type extension
     * is not important for the determination of the file type but it must be specified in the path. Accepted input
     * formats: MDL Molfile, MDL Structure data file (SDF) and SMILES files (of format: [SMILES string][space][name]
     * in each line, see example file). The final test for whether the file is suitable is done in execute(). This option
     * and its argument are always required.
     * <br>* option -t --typeOfMoietiesToRemove <integer>: A number of ["1","2","3"] indicating which type of sugar
     * moieties should be removed, "1" for circular sugar moieties, "2" for linear sugar moieties, or "3" for circular
     * AND linear sugar moieties. Check isLegalTypeOfMoietiesToRemove(int) for the correct mapping of value to option.
     * This option and its argument are always required.
     * <p>
     * * option -glyBond --detectCircSugOnlyWithGlyBond <boolean>: Either "true" or "false", indicating whether circular
     * sugars should be detected (and removed) only if they have an O-glycosidic bond to another moiety or the core of
     * the molecule. Any other value of this argument will be interpreted as "false". Default: "false". This option
     * is optional.
     * <br>* option -remTerm,--removeOnlyTerminalSugars <boolean>: Either "true" or "false", indicating whether only
     * terminal sugar moieties should be removed. Any other value of this argument will be interpreted as "false".
     * Default: "true". Important note: If this setting is set to "true", the input molecules must all consist of one
     * connected structure, respectively. If they already contain multiple, disconnected structures (e.g. counter-ions),
     * the respective molecules are ignored. This option is optional.
     * <br>* option -presMode --preservationModeOption <integer>: A number of ["0","1","2"] indicating which preservation
     * mode to use. This specifies under what circumstances to discard structures that get disconnected from the central
     * core in the sugar removal process, "0" to preserve all disconnected structures (note: this might lead to no
     * circular sugar moieties being detected, depending on the other settings), "1" to remove disconnected structures
     * that do not have enough heavy atoms, or "2" to remove disconnected structures that do not have a sufficient
     * molecular weight. Default: "1" (judge disconnected structures by their heavy atom count). check
     * SugarRemovalUtility.PreservationModeOption enum for the correct mapping of value to option, it corresponds to the
     * ordinal values of the enum constants. This option is optional.
     * <br>* option -presThres --preservationModeThreshold <integer>: An integer number giving the threshold of the
     * preservation mode, i.e. how many heavy atoms a disconnected structure needs to have at least to be not removed or
     * how heavy (in terms of its molecular weight) it needs to be. Default: "5" (heavy atoms). The integer number must
     * be positive. If the option -presMode --preservationModeOption was passed the value "0" (preserve all structures), this
     * option must also be passed a zero value. In the opposite case, this option must be passed a non-zero value
     * if the option -presMode --preservationModeOption was given the value 1 or 2. This option is optional.
     * <br>* option -oxyAtoms --detectCircSugOnlyWithEnoughExocycOxyAtoms <boolean>: Either "true" or "false", indicating
     * whether circular sugars should be detected (and removed) only if they have a sufficient number of attached
     * exocyclic oxygen atoms. Any other value of this argument will be interpreted as "false". Default: "true".
     * This option is optional.
     * <br>* option -oxyAtomsThres --exocycOxyAtomsToAtomsInRingRatioThreshold <number>: A number giving the minimum
     * attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be detected as such.
     * Default: "0.5" (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms).
     * If the option -oxyAtoms --detectCircSugOnlyWithEnoughExocycOxyAtoms was passed the value "false" (detect circular sugars
     * neglecting their number of attached exocyclic oxygen atoms), this option must be passed a zero
     * value. In the opposite case, this option must be passed a non-zero value if the option -oxyAtoms
     * --detectCircSugOnlyWithEnoughExocycOxyAtoms was given the value "true". The number must be positive. This option
     * is optional.
     * <br>* option -linSugInRings --detectLinSugInRings <boolean>: Either "true" or "false", indicating whether linear
     * sugars in rings should be detected (and removed). Any other value of this argument will be interpreted as "false".
     * Default: "false". This option is optional.
     * <br>* option -linSugMinSize --linSugCandidateMinimumSize <integer>: An integer number indicating the minimum
     * number of carbon atoms a linear sugar needs to have to be detected as such. Default: "4". The integer number must
     * be positive and higher than or equal to 1 and also smaller than the linear sugar candidate maximum size
     * (option -linSugMaxSize --linSugCandidateMaximumSize). This option is optional.
     * <br>* option -linSugMaxSize --linSugCandidateMaximumSize <integer>: An integer number indicating the maximum
     * number of carbon atoms a linear sugar needs to have to be detected as such. Default: "7". The integer number must
     * be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size
     * (option -linSugMinSize --linSugCandidateMinimumSize). This option is optional.
     * <br>* option -linAcSug --detectLinAcidicSug <boolean>: Either "true" or "false", indicating whether linear acidic
     * sugars should be included in the set of linear sugar patterns for the initial detection. Any other value of this
     * argument will be interpreted as "false". Default: "false". This option is optional.
     * <br>* option -circSugSpiro --detectSpiroRingsAsCircSug <boolean>: Either "true" or "false", indicating whether
     * spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection.
     * Any other value of this argument will be interpreted as "false". Default: "false". This option is optional.
     * <br>* option -circSugKetoGroups --detectCircularSugarsWithKetoGroups <boolean>: Either "true" or "false",
     * indicating whether circular sugar-like moieties with keto groups should be detected. Any other value of this
     * argument will be interpreted as "false". Default: "false". This option is optional.
     * <p>
     * Example (all settings in default): String[] args = new String[] {"-i", "smiles_test_file.txt", "-t", "3",
     * "-glyBond", "false", "-remTerm", "true", "-presMode", "1", "-presThres", "5", "-oxyAtoms", "true", "-oxyAtomsThres",
     * "0.5", "-linSugInRings", "false", "-linSugMinSize", "4", "-linSugMaxSize", "7", "-linAcSug", "false",
     * "-circSugSpiro", "false", "-circSugKetoGroups", "false"};
     *
     * @throws IllegalArgumentException if any parameter does not meet the specified requirements
     */
    public SugarRemovalUtilityCmdApplication(String[] args) throws IllegalArgumentException {
        CommandLineParser tmpCmdParser = new DefaultParser();
        CommandLine tmpCommandLine;
        try {
            tmpCommandLine = tmpCmdParser.parse(SugarRemovalUtilityCmdApplication.HELP_AND_VERSION_OPTIONS, args, true);
        } catch (ParseException aParseException) {
            throw new IllegalArgumentException(aParseException.getMessage());
        }
        this.wasHelpOrVersionQueried = false;
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.HELP_OPTION.getOpt())) {
            this.wasHelpOrVersionQueried = true;
            HelpFormatter tmpHelpFormatter = new HelpFormatter();
            tmpHelpFormatter.setOptionComparator(null);
            tmpHelpFormatter.setWidth(160);
            System.out.println();
            tmpHelpFormatter.printHelp("java -jar SugarRemovalUtility-jar-with-dependencies.jar",
                    SugarRemovalUtilityCmdApplication.CMD_OPTIONS, true);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.VERSION_OPTION.getOpt())) {
            this.wasHelpOrVersionQueried = true;
            System.out.println();
            System.out.println("Sugar Removal Utility Command Line Application version " + SugarRemovalUtilityCmdApplication.VERSION);
        }
        if (this.wasHelpOrVersionQueried) {
            return;
        }
        tmpCmdParser = new DefaultParser();
        try {
            tmpCommandLine = tmpCmdParser.parse(SugarRemovalUtilityCmdApplication.CMD_OPTIONS, args, false);
        } catch (ParseException aParseException) {
            throw new IllegalArgumentException(aParseException.getMessage());
        }
        if (!tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getOpt())) {
            throw new IllegalArgumentException("Option -" + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getOpt()
                    + " or --" + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getLongOpt()
                    + " is always required but has not been given.");
        }
        String tmpPath = tmpCommandLine.getOptionValue(SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getOpt()).trim();
        if (Objects.isNull(tmpPath) || tmpPath.isBlank()) {
            throw new IllegalArgumentException("Given path to SMILES, SD, or MOL file (option -"
                    + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getOpt() + " or --"
                    + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getLongOpt()
                    + ") is empty or blank.");
        }
        File tmpFile = new File(tmpPath);
        if (!tmpFile.exists() || !tmpFile.canRead() || !tmpFile.isFile()) {
            throw new IllegalArgumentException("Given path to SMILES, SD, or MOL file (option -"
                    + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getOpt() + " or --"
                    + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getLongOpt()
                    + ") does not exist or the file cannot be read.");
        }
        this.inputFile = tmpFile;
        //determination of file type is done in execute(), might cause exceptions there!
        if (!tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION.getOpt())) {
            throw new IllegalArgumentException("Option -" + SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION.getOpt()
                    + " or --" + SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION.getLongOpt()
                    + " is always required but has not been given.");
        }
        int tmpTypeOfMoietiesToRemove = -1;
        try {
            tmpTypeOfMoietiesToRemove = Integer.parseInt(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION.getOpt()).trim());
        } catch (NumberFormatException | NullPointerException anException) {
            throw new IllegalArgumentException("Integer indicating which type of moieties is to remove (option -"
                    + SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION.getOpt() + " or --"
                    + SugarRemovalUtilityCmdApplication.TYPE_OF_MOIETIES_TO_REMOVE_OPTION.getLongOpt()
                    + ") cannot be parsed.");
        }
        if (!SugarRemovalUtilityCmdApplication.isLegalTypeOfMoietiesToRemove(tmpTypeOfMoietiesToRemove)) {
            throw new IllegalArgumentException("Option -" + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getOpt()
                    + " or --" + SugarRemovalUtilityCmdApplication.INPUT_FILE_PATH_OPTION.getLongOpt()
                    + " indicating which type of moieties is to remove must be "
                    + "either 1 (circular sugar moieties), 2 (linear sugar moieties), or 3 (both).");
        }
        this.typeOfMoietiesToRemove = tmpTypeOfMoietiesToRemove;
        this.sugarRemovalUtil = new SugarRemovalUtility(DefaultChemObjectBuilder.getInstance());
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to detect only circular sugars that "
                        + "have an O-glycosidic bond (option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpDetectCircularSugarsOnlyWithOGlycosidicBondSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_O_GLYCOSIDIC_BOND_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(
                    tmpDetectCircularSugarsOnlyWithOGlycosidicBondSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.REMOVE_ONLY_TERMINAL_SUGARS_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.REMOVE_ONLY_TERMINAL_SUGARS_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to remove only terminal sugar moieties (option -"
                        + SugarRemovalUtilityCmdApplication.REMOVE_ONLY_TERMINAL_SUGARS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.REMOVE_ONLY_TERMINAL_SUGARS_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpRemoveOnlyTerminalSugarsSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.REMOVE_ONLY_TERMINAL_SUGARS_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(
                    tmpRemoveOnlyTerminalSugarsSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt())) {
            int tmpPreservationModeSetting = -1;
            try {
                tmpPreservationModeSetting = Integer.parseInt(tmpCommandLine.getOptionValue(
                        SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt()).trim());
            } catch (NumberFormatException | NullPointerException anException) {
                throw new IllegalArgumentException("Integer indicating the preservation mode setting (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getLongOpt()
                        + ") cannot be parsed.");
            }
            //the given int is taken as ordinal value of the indicated enum object; ordinals start at 0
            if (tmpPreservationModeSetting >= SugarRemovalUtility.PreservationModeOption.values().length ||
                    tmpPreservationModeSetting < 0) {
                throw new IllegalArgumentException("Integer indicating the preservation mode setting (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getLongOpt()
                        + ") must be either 0 (keep all), 1 (judge by heavy atom count), or 2 (judge by molecular weight).");
            }
            SugarRemovalUtility.PreservationModeOption tmpPreservationMode = null;
            for (SugarRemovalUtility.PreservationModeOption tmpEnumConstant : SugarRemovalUtility.PreservationModeOption.values()) {
                if (tmpPreservationModeSetting == tmpEnumConstant.ordinal()) {
                    tmpPreservationMode = tmpEnumConstant;
                }
            }
            if (Objects.isNull(tmpPreservationMode)) {
                throw new IllegalArgumentException("Given preservation mode setting (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getLongOpt()
                        + ") does not correspond to an ordinal in the PreservationModeOption enumeration.");
            }
            this.sugarRemovalUtil.setPreservationModeSetting(tmpPreservationMode);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getOpt())) {
            int tmpPreservationModeThresholdSetting = -1;
            try {
                tmpPreservationModeThresholdSetting = Integer.parseInt(tmpCommandLine.getOptionValue(
                        SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getOpt()).trim());
            } catch (NumberFormatException | NullPointerException anException) {
                throw new IllegalArgumentException("Integer indicating the preservation mode setting threshold (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getLongOpt()
                        + ") cannot be parsed.");
            }
            if (tmpPreservationModeThresholdSetting < 0) {
                throw new IllegalArgumentException("Integer indicating the preservation mode setting threshold (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getLongOpt()
                        + ") cannot be negative.");
            }
            //case 1: Preservation mode setting is 'all' (ordinal 0). Then, a 0 threshold should be passed.
            // case 2: Preservation mode is not 'all' (ordinal nonzero). Then, a nonzero threshold should be passed.
            SugarRemovalUtility.PreservationModeOption tmpPreservationModeSetting = this.sugarRemovalUtil.getPreservationModeSetting();
            if (tmpPreservationModeSetting.equals(SugarRemovalUtility.PreservationModeOption.ALL) && tmpPreservationModeThresholdSetting != 0) {
                throw new IllegalArgumentException("Preservation mode 'all' was selected or used as default (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getLongOpt()
                        + "). Therefore, passing a nonzero threshold (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getLongOpt()
                        + ") makes no sense.");
            } else if (!tmpPreservationModeSetting.equals(SugarRemovalUtility.PreservationModeOption.ALL) && tmpPreservationModeThresholdSetting == 0) {
                throw new IllegalArgumentException("Please select preservation mode 'all' (set option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_OPTION.getLongOpt()
                        + " to 1) if the associated threshold should be 0 (option -"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.PRESERVATION_MODE_THRESHOLD_OPTION.getLongOpt()
                        + ").");
            }
            this.sugarRemovalUtil.setPreservationModeThresholdSetting(tmpPreservationModeThresholdSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to detect circular sugars only with " +
                        "a sufficient number of exocyclic oxygen atoms attached (option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(
                    tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getOpt())) {
            double tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting = -1;
            try {
                tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting = Double.parseDouble(tmpCommandLine.getOptionValue(
                        SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getOpt()).trim());
            } catch (NumberFormatException | NullPointerException anException) {
                throw new IllegalArgumentException("Number indicating the exocyclic oxygen atoms to atoms in ring ratio threshold (option -"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getLongOpt()
                        + ") cannot be parsed.");
            }
            boolean tmpIsFinite = Double.isFinite(tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting);
            boolean tmpIsNegative = (tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting < 0);
            if (!tmpIsFinite || tmpIsNegative) {
                throw new IllegalArgumentException("Number indicating the exocyclic oxygen atoms to atoms in ring ratio threshold (option -"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getLongOpt()
                        + ") is NaN, infinite or negative.");
            }
            //case 1: If the number of exocyclic oxygen atoms is neglected, a 0 threshold should be passed
            // case 2: If it is detected, a nonzero threshold should be passed
            boolean tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting =
                    this.sugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected();
            if (!tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
                    && tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting != 0) {
                throw new IllegalArgumentException("The number of exocyclic oxygen atoms of circular sugars is neglected at detection (option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getLongOpt()
                        + " is set to 'false' or 'false' was used as default option). "
                        + "Therefore, passing a nonzero threshold (option -"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getLongOpt()
                        + ") makes no sense.");
            } else if (tmpDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting
                    && tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting == 0) {
                throw new IllegalArgumentException("Please select to neglect the number of exocyclic oxygen atoms (set option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_ONLY_WITH_ENOUGH_EXOCYCLIC_OXYGEN_ATOMS_OPTION.getLongOpt()
                        + " to false) if the associated threshold should be 0 (option -"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.EXOCYCLIC_OXYGEN_ATOMS_TO_ATOMS_IN_RING_RATIO_THRESHOLD_OPTION.getLongOpt()
                        + ").");
            }
            this.sugarRemovalUtil.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(tmpExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.DETECT_LINEAR_SUGARS_IN_RINGS_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_LINEAR_SUGARS_IN_RINGS_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to detect linear sugars inside of rings (option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_LINEAR_SUGARS_IN_RINGS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_LINEAR_SUGARS_IN_RINGS_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpDetectLinearSugarsInRingsSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.DETECT_LINEAR_SUGARS_IN_RINGS_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setDetectLinearSugarsInRingsSetting(
                    tmpDetectLinearSugarsInRingsSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getOpt()) ||
                tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getOpt())) {
            int tmpLinearSugarCandidateMinSizeSetting = -1;
            int tmpLinearSugarCandidateMaxSizeSetting = -1;
            if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getOpt())) {
                try {
                    tmpLinearSugarCandidateMinSizeSetting = Integer.parseInt(tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getOpt()).trim());
                } catch (NumberFormatException | NullPointerException anException) {
                    throw new IllegalArgumentException("The linear sugar candidate minimum size (option -"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getOpt() + " or --"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getLongOpt()
                            + ") cannot be parsed.");
                }
                if (tmpLinearSugarCandidateMinSizeSetting < 1) {
                    throw new IllegalArgumentException("The linear sugar candidate minimum size (option -"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getOpt() + " or --"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getLongOpt()
                            + ") must be equal or higher than 1.");
                }
            }
            if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getOpt())) {
                try {
                    tmpLinearSugarCandidateMaxSizeSetting = Integer.parseInt(tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getOpt()).trim());
                } catch (NumberFormatException | NullPointerException anException) {
                    throw new IllegalArgumentException("The linear sugar candidate maximum size (option -"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getOpt() + " or --"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getLongOpt()
                            + ") cannot be parsed.");
                }
                if (tmpLinearSugarCandidateMaxSizeSetting < 1) {
                    throw new IllegalArgumentException("The linear sugar candidate maximum size (option -"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getOpt() + " or --"
                            + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getLongOpt()
                            + ") must be equal or higher than 1.");
                }
            }
            if (tmpLinearSugarCandidateMinSizeSetting == -1) {
                tmpLinearSugarCandidateMinSizeSetting = this.sugarRemovalUtil.getLinearSugarCandidateMinSizeSetting();
            }
            if (tmpLinearSugarCandidateMaxSizeSetting == -1) {
                tmpLinearSugarCandidateMaxSizeSetting = this.sugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting();
            }
            if (tmpLinearSugarCandidateMinSizeSetting > tmpLinearSugarCandidateMaxSizeSetting) {
                throw new IllegalArgumentException("The linear sugar candidate minimum size (option -"
                        + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_OPTION.getLongOpt()
                        + ", set to " + tmpLinearSugarCandidateMinSizeSetting + ", default value: "
                        + SugarRemovalUtility.LINEAR_SUGAR_CANDIDATE_MIN_SIZE_DEFAULT
                        + ") must be smaller than the maximum size (option -"
                        + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_OPTION.getLongOpt()
                        + ", set to " + tmpLinearSugarCandidateMaxSizeSetting + ", default value: "
                        + SugarRemovalUtility.LINEAR_SUGAR_CANDIDATE_MAX_SIZE_DEFAULT
                        + ").");
            }
            this.sugarRemovalUtil.setLinearSugarCandidateMinSizeSetting(tmpLinearSugarCandidateMinSizeSetting);
            this.sugarRemovalUtil.setLinearSugarCandidateMaxSizeSetting(tmpLinearSugarCandidateMaxSizeSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.DETECT_LINEAR_ACIDIC_SUGARS_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_LINEAR_ACIDIC_SUGARS_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to include linear acidic sugars in " +
                        "the patterns for initial detection of linear sugars(option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_LINEAR_ACIDIC_SUGARS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_LINEAR_ACIDIC_SUGARS_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpDetectLinearAcidicSugarsSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.DETECT_LINEAR_ACIDIC_SUGARS_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setDetectLinearAcidicSugarsSetting(
                    tmpDetectLinearAcidicSugarsSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to detect spiro rings as circular sugars (option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpDetectSpiroRingsAsCircularSugarsSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.DETECT_SPIRO_RINGS_AS_CIRCULAR_SUGARS_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(
                    tmpDetectSpiroRingsAsCircularSugarsSetting);
        }
        if (tmpCommandLine.hasOption(SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION.getOpt())) {
            if (Objects.isNull(tmpCommandLine.getOptionValue(
                    SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION.getOpt()))) {
                throw new IllegalArgumentException("Boolean value indicating whether to detect circular sugars with keto groups (option -"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION.getOpt() + " or --"
                        + SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION.getLongOpt()
                        + ") is empty (null).");
            }
            boolean tmpDetectCircularSugarsWithKetoGroupsSetting = Boolean.parseBoolean(
                    tmpCommandLine.getOptionValue(
                            SugarRemovalUtilityCmdApplication.DETECT_CIRCULAR_SUGARS_WITH_KETO_GROUPS_OPTION.getOpt()).trim());
            this.sugarRemovalUtil.setDetectCircularSugarsWithKetoGroupsSetting(
                    tmpDetectCircularSugarsWithKetoGroupsSetting);
        }
        this.sugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Returns true if the -h --help or -v --version options were used at the command-line. In this case, the object
     * is not properly instantiated and the application should be exited.
     *
     * @return true if the application should be exited
     */
    public boolean wasHelpOrVersionQueried() {
        return this.wasHelpOrVersionQueried;
    }

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
        if (Objects.isNull(this.inputFile) || Objects.isNull(this.sugarRemovalUtil) || this.typeOfMoietiesToRemove == 0) {
            throw new IllegalArgumentException("This object has not been properly instantiated because the usage " +
                    "instructions or version were queried by the command-line. It therefore cannot be executed.");
        }
        System.out.println();
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
        System.out.println();
        System.out.println("Given settings for Sugar Removal Utility:");
        switch (this.typeOfMoietiesToRemove) {
            case 1:
                System.out.println("\tType of moieties to remove: circular sugars (1)");
                break;
            case 2:
                System.out.println("\tType of moieties to remove: linear sugars (2)");
                break;
            case 3:
                System.out.println("\tType of moieties to remove: circular AND linear sugars (3)");
                break;
            default:
                throw new IllegalArgumentException("Type of moieties to remove must be either 1 (circular sugar moieties), " +
                        "2 (linear sugar moieties), or 3 (both).");
        }
        System.out.println("\tDetect circular sugars only with O-glycosidic bond: "
                + this.sugarRemovalUtil.areOnlyCircularSugarsWithOGlycosidicBondDetected());
        System.out.println("\tRemove only terminal sugar moieties: " + this.sugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println("\tPreservation mode setting: " + this.sugarRemovalUtil.getPreservationModeSetting().name()
                + " (" + this.sugarRemovalUtil.getPreservationModeSetting().ordinal() + ")");
        System.out.println("\tPreservation mode threshold: "
                + this.sugarRemovalUtil.getPreservationModeThresholdSetting());
        System.out.println("\tDetect circular sugars only with enough exocyclic oxygen atoms: "
                + this.sugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected());
        System.out.println("\tExocyclic oxygen atoms to atoms in ring ratio threshold: "
                + this.sugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        System.out.println("\tDetect linear sugars in rings: " + this.sugarRemovalUtil.areLinearSugarsInRingsDetected());
        System.out.println("\tLinear sugar candidate minimum size: "
                + this.sugarRemovalUtil.getLinearSugarCandidateMinSizeSetting() + " carbon atoms");
        System.out.println("\tLinear sugar candidate maximum size: "
                + this.sugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting() + " carbon atoms");
        System.out.println("\tDetect linear acidic sugars: " + this.sugarRemovalUtil.areLinearAcidicSugarsDetected());
        System.out.println("\tDetect spiro rings as circular sugars: "
                + this.sugarRemovalUtil.areSpiroRingsDetectedAsCircularSugars());
        System.out.println("\tDetect circular sugars with keto groups: "
                + this.sugarRemovalUtil.areCircularSugarsWithKetoGroupsDetected());
        System.out.println();
        System.out.println("Entering iteration of molecules...");
        System.out.println();
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
        System.out.println();
        System.out.println((tmpMoleculeCounter - 1) + " molecules were analysed.");
        System.out.println(tmpFatalExceptionCounter + " fatal exceptions occurred (molecule could not be read, was " +
                "disconnected or could not be cloned).");
        System.out.println(tmpMinorExceptionsCounter + " minor exceptions occurred (SMILES codes could not be generated etc.).");
        System.out.println(tmpSugarContainingMoleculesCounter + " molecules contained sugar moieties.");
        System.out.println(tmpContainsCircularSugarsCounter + " molecules contained circular sugar moieties.");
        System.out.println(tmpContainsLinearSugarsCounter + " molecules contained linear sugar moieties.");
        System.out.println(tmpBasicallyASugarCounter + " molecules were basically sugars because they were completely removed.");
        System.out.println();
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
     * input value for this option.
     *
     * @param aTypeOfMoietiesToRemove an integer that is supposed to be mapped to a type of moiety to remove
     * @return true if the given value is valid for the respective option
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
