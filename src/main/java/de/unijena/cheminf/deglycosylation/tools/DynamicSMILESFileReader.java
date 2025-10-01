/*
 * MIT License
 *
 * Copyright (c) 2025 Jonas Schaub, Achim Zielesny, Christoph Steinbeck, Maria Sorokina
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

package de.unijena.cheminf.deglycosylation.tools;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.io.iterator.DefaultIteratingChemObjectReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * File reader for different kinds of files with a SMILES code column. The reader can detect the structure of the file based
 * on the first few lines and some assumptions like the SMILES and ID/name columns should be the first two but can
 * be in both positions. Unsuitable for reaction SMILES or CxSMILES.
 *
 * @author Jonas Schaub
 * @author Samuel Behr
 * @version 1.0.0.0
 */
public class DynamicSMILESFileReader extends DefaultIteratingChemObjectReader<IAtomContainer> {
    //<editor-fold desc="Public static final class constants">
    /**
     * Possible SMILES file separators used to separate SMILES code from ID. Ordered so that non-whitespace characters
     * are tested first.
     */
    public static final Set<String> POSSIBLE_SMILES_FILE_SEPARATORS = Set.of(",", ";", " ", "\t");
    //
    /**
     * Maximum number of lines starting from the first one to check for valid SMILES strings in a SMILES file when
     * trying to determine the SMILES code column and separator.
     */
    public static final int MAXIMUM_LINE_NUMBER_TO_CHECK_IN_SMILES_FILES = 10;
    //
    /**
     * Strings that can be parsed by CDK SmilesParser as SMILES codes but should be ignored when detecting the file structure,
     * e.g. "ID" is a likely column name but could be parsed into Iodine-Deuterium as a SMILES code.
     */
    public static final Set<String> PARSABLE_SMILES_EXCEPTIONS = Set.of("ID");
    //</editor-fold>
    //
    //<editor-fold desc="Private static class constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(DynamicSMILESFileReader.class.getName());
    //
    /**
     * Parser for SMILES codes.
     */
    private final SmilesParser smilesParser;
    //</editor-fold>
    //
    //<editor-fold desc="Private variables">
    /**
     * Number of lines that were skipped (empty or erroneous SMILES codes etc.) in the last file import of this instance,
     * headline not counted.
     */
    private int skippedLinesCounter;
    //
    /**
     * Buffered reader to read the file line by line.
     */
    private BufferedReader buffReader;
    //
    /**
     * Format of the SMILES file to read.
     */
    private DynamicSMILESFileFormat format;
    //
    /**
     * Current line read from the file, null if the next line has not been read yet.
     */
    private String currentLine;
    //
    /**
     * Counter for the current line in the file, starting at -1 before reading the first line (0 after reading the
     * headline if present).
     */
    private int lineInFileCounter;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructs a new DynamicSMILESFileReader that can read molecules from a given a
     * InputStream.
     *
     * @param in the {@link InputStream} to read from
     */
    public DynamicSMILESFileReader(InputStream in, DynamicSMILESFileFormat aFormat) throws CDKException {
        this(new InputStreamReader(in), aFormat);
    }
    //
    /**
     * Constructs a new DynamicSMILESFileReader that can read molecules from a given file.
     */
    public DynamicSMILESFileReader(File aFile, DynamicSMILESFileFormat aFormat) throws CDKException, FileNotFoundException {
        this(new FileReader(aFile), aFormat);
    }
    //
    /**
     * Constructs a new DynamicSMILESFileReader that can read molecules from a given a
     * Reader.
     *
     * @param in the {@link Reader} to read from
     */
    public DynamicSMILESFileReader(Reader in, DynamicSMILESFileFormat aFormat) throws CDKException {
        this.setSMILESFileFormat(aFormat);
        this.setReader(in);
        this.smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the number of lines that were skipped (empty or erroneous SMILES codes etc.) in the last file import of
     * this instance (reset at each new import), headline not counted.
     *
     * @return nr of lines skipped in last import
     */
    public int getSkippedLinesCounter() {
        return skippedLinesCounter;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    /**
     * Checking the first few lines of a SMILES file for parsable SMILES codes and saving the determined separator
     * character and SMILES code and ID column positions.
     * Expects one parsable SMILES code per line of the file and an optional second element, which is interpreted as the
     * molecule's ID or name and is separated from the SMILES code by one of the separator tokens tab, semicolon, comma, or space.
     * Unsuitable for reaction SMILES or CxSMILES.
     *
     * @param lines first few lines of a SMILES file
     * @return determined format of the given lines
     * @throws IOException if the lines do not adhere to the format expectations
     */
    public static DynamicSMILESFileFormat detectFormat(List<String> lines) throws IOException {
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        // AtomContainer to save the parsed SMILES in
        IAtomContainer tmpMolecule = tmpBuilder.newAtomContainer();
        SmilesParser tmpSmilesParser = new SmilesParser(tmpBuilder);
        String tmpSmilesFileDeterminedSeparator = String.valueOf(DynamicSMILESFileFormat.PLACEHOLDER_SEPARATOR_CHAR);
        int tmpSmilesCodeExpectedPosition = DynamicSMILESFileFormat.DEFAULT_SMILES_COLUMN_POSITION;
        int tmpIDExpectedPosition = DynamicSMILESFileFormat.PLACEHOLDER_ID_COLUMN_POSITION;
        int tmpCurrentLineInFileCounter = -1;
        String tmpSmilesFileFirstLine = "";
        findSeparatorLoop:
        for (String tmpSmilesFileCurrentLine : lines) {
            tmpCurrentLineInFileCounter++;
            if (Objects.isNull(tmpSmilesFileCurrentLine)) {
                break findSeparatorLoop;
            }
            if (tmpCurrentLineInFileCounter == 0) {
                // saved for determination of whether the file has a headline below
                tmpSmilesFileFirstLine = tmpSmilesFileCurrentLine;
            }
            // first try the whole line because the file might have only one column
            if (!tmpSmilesFileCurrentLine.trim().isEmpty()
                    //not trimmed because whitespaces are invalid and should be detected
                    && DynamicSMILESFileReader.containsOnlySMILESValidCharacters(tmpSmilesFileCurrentLine)
                    && !DynamicSMILESFileReader.PARSABLE_SMILES_EXCEPTIONS.contains(tmpSmilesFileCurrentLine.trim())) {
                try {
                    //if parsing fails goes to catch block below
                    // trimmed because a leading whitespace followed by a character string are interpreted as an empty structure and its name
                    tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmilesFileCurrentLine.trim());
                    if (!tmpMolecule.isEmpty()) {
                        //success, SMILES column is identified
                        tmpSmilesCodeExpectedPosition = 0;
                        break findSeparatorLoop;
                    }
                } catch (InvalidSmilesException anException) {
                    // do nothing, continue with splitting the line using different separators
                }
            }
            for (String tmpSeparator : DynamicSMILESFileReader.POSSIBLE_SMILES_FILE_SEPARATORS) {
                // limit param = 3 because we assume that SMILES code and ID are  in the first two columns
                String[] tmpProcessedLineArray = tmpSmilesFileCurrentLine.split(tmpSeparator, 3);
                for (int i = 0; i < tmpProcessedLineArray.length; i++) {
                    String tmpNextElementOfLine = tmpProcessedLineArray[i];
                    if (tmpNextElementOfLine.trim().isEmpty()
                            || !DynamicSMILESFileReader.containsOnlySMILESValidCharacters(tmpNextElementOfLine)
                            || DynamicSMILESFileReader.PARSABLE_SMILES_EXCEPTIONS.contains(tmpNextElementOfLine.trim())) {
                        continue; //... to try next element in row
                    }
                    // check only the first two columns, see comment above
                    if (i > 1) {
                        break; //... try next separator
                    }
                    try {
                        //if parsing fails goes to catch block below
                        tmpMolecule = tmpSmilesParser.parseSmiles(tmpNextElementOfLine.trim());
                        if (!tmpMolecule.isEmpty()) {
                            //success, separator and SMILES column are identified
                            tmpSmilesFileDeterminedSeparator = tmpSeparator;
                            tmpSmilesCodeExpectedPosition = i;
                            if (tmpProcessedLineArray.length > 1) {
                                tmpIDExpectedPosition = tmpSmilesCodeExpectedPosition == 0 ? 1 : 0;
                            }
                            break findSeparatorLoop;
                        } //else {
                        // continue to try next element in row
                        //}
                    } catch (InvalidSmilesException anException) {
                        // continue to try next element in row
                    }
                }
            }
        }
        if (tmpMolecule.isEmpty()) {
            throw new IOException("Chosen file does not fit to the expected format of a SMILES file.");
        }
        boolean tmpHasHeaderLine;
        try {
            if (tmpIDExpectedPosition == -1) {
                tmpSmilesParser.parseSmiles(tmpSmilesFileFirstLine.trim());
            } else {
                tmpSmilesParser.parseSmiles(tmpSmilesFileFirstLine.trim().split(tmpSmilesFileDeterminedSeparator, 3)[tmpSmilesCodeExpectedPosition]);
            }
            tmpHasHeaderLine = false;
        } catch (InvalidSmilesException anException) {
            tmpHasHeaderLine = true;
        }
        return new DynamicSMILESFileFormat(tmpHasHeaderLine, tmpSmilesFileDeterminedSeparator.charAt(0), tmpSmilesCodeExpectedPosition, tmpIDExpectedPosition);
    }
    //
    /**
     * Checking the first few lines of the SMILES file for parsable SMILES codes and saving the determined separator
     * character and SMILES code and ID column positions.
     * Expects one parsable SMILES code per line of the file and an optional second element, which is interpreted as the
     * molecule's ID or name and is separated from the SMILES code by one of the separator tokens tab, semicolon, comma, or space.
     * Unsuitable for reaction SMILES or CxSMILES.
     *
     * @param aFile a SMILES file
     * @return determined format of the given file
     * @throws IOException if the file cannot be found or does not adhere to the format expectations
     */
    public static DynamicSMILESFileFormat detectFormat(File aFile) throws IOException {
        try (
                // throws FileNotFoundException if file cannot be found, see catch block below
                FileReader tmpSmilesFileReader = new FileReader(aFile);
                BufferedReader tmpSmilesFileBufferedReader = new BufferedReader(tmpSmilesFileReader, 65536)
        ) {
            String tmpCurrentLine;
            int tmpLinesCheckedCounter = 0;
            List<String> tmpLinesToCheck = new ArrayList<>(DynamicSMILESFileReader.MAXIMUM_LINE_NUMBER_TO_CHECK_IN_SMILES_FILES);
            while (tmpLinesCheckedCounter < DynamicSMILESFileReader.MAXIMUM_LINE_NUMBER_TO_CHECK_IN_SMILES_FILES
                    && (tmpCurrentLine = tmpSmilesFileBufferedReader.readLine()) != null) {
                tmpLinesToCheck.add(tmpCurrentLine);
                tmpLinesCheckedCounter++;
            }
            return DynamicSMILESFileReader.detectFormat(tmpLinesToCheck);
        } catch (FileNotFoundException anException) {
            String tmpMessage = "File " + aFile.getPath() + " could not be found";
            DynamicSMILESFileReader.LOGGER.log(Level.SEVERE, tmpMessage);
            throw new IOException(tmpMessage);
        }
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public methods">
    /**
     * Reads SMILES file according to the given format. Splits the lines at the given separator character, ignores the
     * first line if the format defines that the file has a headline, parses SMILES codes and IDs from the defined columns, etc.
     * Skipped lines (due to being empty or containing erroneous SMILES codes) are counted and this counter can be queried
     * after import via the respective getter method. If a name/ID column is given in the file, it is read and saved as
     * a property of the respective atom container under the name property key taken from the Importer class.
     *
     * @return atom container set parsed from the file
     * @throws IOException if the given file cannot be found
     */
    public IAtomContainerSet readToSet() throws IOException {
        IAtomContainerSet tmpAtomContainerSet = new AtomContainerSet();
        while (this.hasNext()) {
            IAtomContainer next = this.next();
            if (next != null)
                tmpAtomContainerSet.addAtomContainer(next);
        }
        return tmpAtomContainerSet;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Package private methods">
    /**
     * Check the given String for characters that are not defined in SMILES encoding. The allowed characters are
     * 0-9 (rings, hydrogen counts, charge counts, or isotopes),
     * a-z, A-Z (element symbols),
     * * (wildcard atoms),
     * [, ] (inorganic atoms or explicit environments),
     * -, + (charges or - for explicit single bonds),
     * @ (tetrahedral stereochemistry),
     * =, #, $ (bonds up to quadruple),
     * % (multi-digit ring numbers),
     * : (tautomer bond),
     * . (disconnected parts),
     * (, ) (branches),
     * /, \ (cis/trans stereochemistry).
     * All other characters, including whitespace characters, are not allowed.
     *
     * @param aPotentialSMILESString the string to test
     * @return true if the input string contains only characters defined in the SMILES format
     */
    static boolean containsOnlySMILESValidCharacters(String aPotentialSMILESString) {
        Pattern pattern = Pattern.compile("^[0-9a-zA-Z*\\[\\]\\-+@=#$%:.()/\\\\]+$");
        Matcher matcher = pattern.matcher(aPotentialSMILESString);
        return matcher.find();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Overwritten DefaultIteratingChemObjectReader methods">

    @Override
    public boolean hasNext() {
        if (this.currentLine == null) {
            try {
                this.currentLine = this.buffReader.readLine();
                return this.currentLine != null;
            } catch (IOException e) {
                return false;
            }
        } else {
            return true;
        }
    }

    @Override
    public IAtomContainer next() {
        if (this.currentLine == null && !this.hasNext()) {
            return null;
        }
        IAtomContainer tmpMolecule;
        String[] tmpProcessedLineArray = new String[0];
        this.lineInFileCounter++;
        //trying to parse as SMILES code
        try {
            String tmpSmiles;
            if (this.format.hasIDColumn()) {
                tmpProcessedLineArray = this.currentLine.split(this.format.getSeparatorChar().toString(), 3);
                tmpSmiles = tmpProcessedLineArray[this.format.getSMILESCodeColumnPosition()].trim().isBlank() ? null :
                        tmpProcessedLineArray[this.format.getSMILESCodeColumnPosition()].trim();
            } else {
                tmpSmiles = this.currentLine.trim();
            }
            this.currentLine = null;
            if (tmpSmiles != null && !tmpSmiles.isEmpty()) {
                //throws exception if SMILES string is null, goes to catch block
                tmpMolecule = this.smilesParser.parseSmiles(tmpSmiles);
            } else {
                throw new InvalidSmilesException("String is empty");
            }
        } catch (InvalidSmilesException | IndexOutOfBoundsException | NullPointerException anException) {
            this.skippedLinesCounter++;
            DynamicSMILESFileReader.LOGGER.log(Level.WARNING, String.format("Import failed for structure in line (starting at 0):\t%s", this.lineInFileCounter), anException);
            return null;
        }
        //setting the name of the atom container
        String tmpName;
        if (this.format.hasIDColumn() && tmpProcessedLineArray.length > 1 && !tmpProcessedLineArray[this.format.getIDColumnPosition()].trim().isEmpty()) {
            tmpName = tmpProcessedLineArray[this.format.getIDColumnPosition()].trim();
            tmpMolecule.setProperty(CDKConstants.TITLE, tmpName);
        }
        return tmpMolecule;
    }

    @Override
    public void setReader(Reader reader) throws CDKException{
        if (reader instanceof BufferedReader) {
            this.buffReader = (BufferedReader) reader;
        } else {
            this.buffReader = new BufferedReader(reader, 65536);
        }
        this.skippedLinesCounter = 0;
        this.lineInFileCounter = -1;
        if (this.format.hasHeaderLine()) {
            try {
                this.buffReader.readLine();
                this.lineInFileCounter++;
            } catch (IOException anException) {
                throw new CDKException(anException.getMessage());
            }
        }
    }

    @Override
    public void setReader(InputStream inputStream) throws CDKException {
        this.setReader(new InputStreamReader(inputStream));
    }

    public void setSMILESFileFormat(DynamicSMILESFileFormat aFormat) {
        this.format = aFormat;
    }

    @Override
    public IResourceFormat getFormat() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() throws IOException {
        this.buffReader.close();
    }
    //</editor-fold>
}
