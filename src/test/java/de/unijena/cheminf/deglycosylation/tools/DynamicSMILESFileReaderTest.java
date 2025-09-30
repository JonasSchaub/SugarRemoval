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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.File;
/**
 * Test class for the DynamicSMILESFileReader class.
 *
 * @author Samuel Behr
 * @author Jonas Schaub
 * @version 2.0.0.0
 */
class DynamicSMILESFileReaderTest {
    /**
     * Test containsOnlySMILESValidCharacters() for false-positives, e.g. two tab-separated strings, some of which can
     * be interpreted by the CDK SmilesParser (it only parses the first part up to the first whitespace character and
     * does not throw an error but interprets the rest as title of the structure).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void containsOnlySMILESValidCharactersTest() throws Exception {
        Assertions.assertFalse(DynamicSMILESFileReader.containsOnlySMILESValidCharacters("CCCCOCCC\tlfdsklhfdfvdbgvb"));
        Assertions.assertFalse(DynamicSMILESFileReader.containsOnlySMILESValidCharacters("CCCCOCCC lfdsklhfdfvdbgvb"));
        Assertions.assertFalse(DynamicSMILESFileReader.containsOnlySMILESValidCharacters(""));
        Assertions.assertFalse(DynamicSMILESFileReader.containsOnlySMILESValidCharacters("\t"));
        for (String tmpSeparator : DynamicSMILESFileReader.POSSIBLE_SMILES_FILE_SEPARATORS) {
            Assertions.assertFalse(DynamicSMILESFileReader.containsOnlySMILESValidCharacters(tmpSeparator), "was true for " + tmpSeparator);
        }
    }
    //
    /**
     * Test file's specifications:
     * - .txt file
     * - with headline
     * - SMILES code column only (no ID or name)
     * - including some blank lines between
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnOneColumnFileWithBlankLinesAndHeadlineTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileOne.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertTrue(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(0, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertFalse(tmpFormat.hasIDColumn());
        Assertions.assertEquals(DynamicSMILESFileFormat.PLACEHOLDER_SEPARATOR_CHAR, tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - .txt file
     * - with headline
     * - SMILES code column only (no ID or name)
     * - including some blank lines between
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnOneColumnFileWithBlankLinesAndHeadlineTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileOne.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        Assertions.assertEquals(3, tmpMolSet.getAtomContainerCount());
        Assertions.assertEquals(2, tmpReader.getSkippedLinesCounter());
    }
    //
    /**
     * Test file's specifications:
     * - .smi file
     * - no headline
     * - ID first in line
     * - used separator: "\t"
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnTwoColumnFileTabSeparatedAndNoHeadlineTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileTwo.smi").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertFalse(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(0, tmpFormat.getIDColumnPosition());
        Assertions.assertEquals(1, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertTrue(tmpFormat.hasIDColumn());
        Assertions.assertEquals('\t', tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - .smi file
     * - no headline
     * - ID first in line
     * - used separator: "\t"
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnTwoColumnFileTabSeparatedAndNoHeadlineTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileTwo.smi").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        Assertions.assertEquals(5, tmpMolSet.getAtomContainerCount());
        Assertions.assertEquals("CNP0337481", tmpMolSet.getAtomContainer(4).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals(0, tmpReader.getSkippedLinesCounter());
    }
    //
    /**
     * Test file's specifications:
     * - Headline
     * - "NAME" second in line and containing spaces
     * - used separator: ";"
     * - two lines with invalid SMILES code
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnTwoColumnFileSemicolonSeparatedWithHeadlineTwoInvalidLinesTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileThree.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertTrue(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(1, tmpFormat.getIDColumnPosition());
        Assertions.assertEquals(0, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertTrue(tmpFormat.hasIDColumn());
        Assertions.assertEquals(';', tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - Headline
     * - "NAME" second in line and containing spaces
     * - used separator: ";"
     * - two lines with invalid SMILES code
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnTwoColumnFileSemicolonSeparatedWithHeadlineTwoInvalidLinesTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileThree.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        Assertions.assertEquals(3, tmpMolSet.getAtomContainerCount());
        Assertions.assertEquals("Istanbulin A", tmpMolSet.getAtomContainer(1).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals("Valdiazen", tmpMolSet.getAtomContainer(2).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals(2, tmpReader.getSkippedLinesCounter());
    }
    //
    /**
     * Test file's specifications:
     * - one single line only
     * - ID first in line
     * - used separator: " "
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnTwoColumnFileSpaceSeparatedWithOnlyOneLineTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileFour.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertFalse(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(0, tmpFormat.getIDColumnPosition());
        Assertions.assertEquals(1, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertTrue(tmpFormat.hasIDColumn());
        Assertions.assertEquals(' ', tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - one single line only
     * - ID first in line
     * - used separator: " "
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnTwoColumnFileSpaceSeparatedWithOnlyOneLineTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileFour.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        Assertions.assertEquals(1, tmpMolSet.getAtomContainerCount());
        Assertions.assertEquals("CNP0356547", tmpMolSet.getAtomContainer(0).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals(0, tmpReader.getSkippedLinesCounter());
    }
    //
    /**
     * Test file's specifications:
     * - headline and blank line first
     * - three elements per line
     * - SMILES first in line, ID second and a neglectable third element
     * - third element in line
     * - used separator: "\t"
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnThreeColumnFileWithHeadlineTabSeparatedTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileFive.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertTrue(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(1, tmpFormat.getIDColumnPosition());
        Assertions.assertEquals(0, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertTrue(tmpFormat.hasIDColumn());
        Assertions.assertEquals('\t', tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - headline and blank line first
     * - three elements per line
     * - SMILES first in line, ID second and a neglectable third element
     * - third element in line
     * - used separator: "\t"
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnThreeColumnFileWithHeadlineTabSeparatedTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileFive.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        String[] tmpTestFileFiveSmiles = new String[] {"OC=1C=C(O)C=C(C1)C=2OC=3C=CC=CC3C2", "OC=1C=C(O)C(=C(C1)C(C)C(O)C)C"};
        String[] tmpTestFileFiveIDs = new String[] {"CNP0192622", "CNP0262448"};
        int i = 0;
        for (IAtomContainer tmpAtomContainer : tmpMolSet.atomContainers()) {
            Assertions.assertEquals(tmpTestFileFiveSmiles[i],DynamicSMILESFileReaderTest.createUniqueSmiles(tmpAtomContainer, false));
            Assertions.assertEquals(tmpTestFileFiveIDs[i],tmpAtomContainer.getProperty(CDKConstants.TITLE));
            i++;
        }
        Assertions.assertEquals(1, tmpReader.getSkippedLinesCounter());
    }
    //
    /**
     * Test file's specifications:
     * - 51 lines, 50 with structures, 1 header line
     * - ID second in line
     * - used separator: " "
     * - multiple garbage columns after the first 2
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnCOCONUTFileWithMultipleColumnsTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileSix.smi").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertTrue(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(1, tmpFormat.getIDColumnPosition());
        Assertions.assertEquals(0, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertTrue(tmpFormat.hasIDColumn());
        Assertions.assertEquals(' ', tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - 51 lines, 50 with structures, 1 header line
     * - ID second in line
     * - used separator: " "
     * - multiple garbage columns after the first 2
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnCOCONUTFileWithMultipleColumnsTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileSix.smi").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        Assertions.assertEquals(50, tmpMolSet.getAtomContainerCount());
        Assertions.assertEquals("CNP0000001", tmpMolSet.getAtomContainer(0).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals(0, tmpReader.getSkippedLinesCounter());
    }
    //
    /**
     * Test file's specifications:
     * - 37 lines
     * - no headline
     * - no faulty lines
     * - ID first in line, SMILES second
     * - used separator: tab
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFormatDetectionOnChEBIFileTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileSeven.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        Assertions.assertTrue(tmpFormat.hasIDColumn());
        Assertions.assertFalse(tmpFormat.hasHeaderLine());
        Assertions.assertEquals(1, tmpFormat.getSMILESCodeColumnPosition());
        Assertions.assertEquals(0, tmpFormat.getIDColumnPosition());
        Assertions.assertEquals('\t', tmpFormat.getSeparatorChar());
    }
    //
    /**
     * Test file's specifications:
     * - 37 lines
     * - no headline
     * - no faulty lines
     * - ID first in line, SMILES second
     * - used separator: tab
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void smilesFileImportOnChEBIFileTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpResourceFile = new File(tmpClassLoader.getResource("SMILESTestFileSeven.txt").getFile());
        DynamicSMILESFileFormat tmpFormat = DynamicSMILESFileReader.detectFormat(tmpResourceFile);
        DynamicSMILESFileReader tmpReader = new DynamicSMILESFileReader(tmpResourceFile, tmpFormat);
        IAtomContainerSet tmpMolSet = tmpReader.readToSet();
        Assertions.assertEquals(37, tmpMolSet.getAtomContainerCount());
        Assertions.assertEquals("cmnpd_id_10213", tmpMolSet.getAtomContainer(0).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals("cmnpd_id_11687", tmpMolSet.getAtomContainer(36).getProperty(CDKConstants.TITLE));
        Assertions.assertEquals(0, tmpReader.getSkippedLinesCounter());
    }

    /**
     * Creates a unique SMILES string out of the given atom container or returns null, if the creation was not possible.
     * If the SMILES could not be created in the first place, it is retried with a kekulized clone of the given atom
     * container. Aromaticity information is encoded in the returned SMILES string, if there is any given. Unique SMILES
     * codes do NOT encode stereochemistry by default! This can be turned on with the second parameter.
     *
     * @param anAtomContainer atom container the unique SMILES should be created of
     * @param isStereoChemEncoded whether stereochemistry should be encoded
     * @return unique SMILES of the given atom container or 'null' if no creation was possible
     */
    private static String createUniqueSmiles(IAtomContainer anAtomContainer, boolean isStereoChemEncoded) {
        return DynamicSMILESFileReaderTest.createUniqueSmiles(anAtomContainer, isStereoChemEncoded, true);
    }

    /**
     * Creates a unique SMILES string out of the given atom container or returns null, if the creation was not possible.
     * If the SMILES could not be created in the first place, it is retried with a kekulized clone of the given atom
     * container. Unique SMILES codes do NOT encode stereochemistry or aromaticity by default! This can be turned on
     * with the parameters.
     *
     * @param anAtomContainer atom container the unique SMILES should be created of
     * @param isStereoChemEncoded whether stereochemistry should be encoded
     * @param isAromaticityEncoded whether aromaticity should be encoded
     * @return unique SMILES of the given atom container or 'null' if no creation was possible
     */
    private static String createUniqueSmiles(IAtomContainer anAtomContainer, boolean isStereoChemEncoded, boolean isAromaticityEncoded) {
        int tmpFlavor = SmiFlavor.Unique;
        if (isAromaticityEncoded) {
            tmpFlavor = tmpFlavor | SmiFlavor.UseAromaticSymbols;
        }
        if (isStereoChemEncoded && anAtomContainer.stereoElements().iterator().hasNext()) {
            tmpFlavor = tmpFlavor | SmiFlavor.Stereo;
        }
        String tmpSmiles = null;
        try {
            try {
                tmpSmiles = SmilesGenerator.create(anAtomContainer, tmpFlavor, new int[anAtomContainer.getAtomCount()]);
            } catch (CDKException anException) {
                IAtomContainer tmpAtomContainer = anAtomContainer.clone();
                Kekulization.kekulize(tmpAtomContainer);
                tmpSmiles = SmilesGenerator.create(tmpAtomContainer, tmpFlavor, new int[anAtomContainer.getAtomCount()]);
                System.out.println(String.format("Kekulized molecule %s", anAtomContainer.getProperty(CDKConstants.TITLE)));
            }
        } catch (CDKException | NullPointerException | IllegalArgumentException | CloneNotSupportedException | ArrayIndexOutOfBoundsException anException){
            System.out.println(String.format("%s; molecule name: %s", anException.toString(), anAtomContainer.getProperty(CDKConstants.TITLE)));
        }
        return tmpSmiles;
    }
}
