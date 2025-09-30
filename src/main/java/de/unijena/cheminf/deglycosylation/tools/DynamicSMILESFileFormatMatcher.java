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

import org.openscience.cdk.io.formats.AbstractResourceFormat;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.formats.IChemFormatMatcher;
import org.openscience.cdk.io.formats.IResourceFormat;
import org.openscience.cdk.io.formats.SMILESFormat;
import org.openscience.cdk.tools.DataFeatures;

import java.io.IOException;
import java.util.List;

/**
 * TODO
 */
public class DynamicSMILESFileFormatMatcher
        extends AbstractResourceFormat
        implements IChemFormatMatcher {

    private static IResourceFormat instance = null;

    public DynamicSMILESFileFormatMatcher() {}

    public static IResourceFormat getInstance() {
        if (instance == null)
            instance = new DynamicSMILESFileFormatMatcher();
        return instance;
    }

    @Override
    public MatchResult matches(List<String> list) {
        try {
            DynamicSMILESFileFormat format = DynamicSMILESFileReader.detectFormat(list);
            return new MatchResult(true,
                    (IChemFormat)SMILESFormat.getInstance(),
                    Integer.valueOf(1));
        } catch (IOException e) {
            return new MatchResult(false, null, Integer.MAX_VALUE);
        }
    }

    @Override
    public String getReaderClassName() {
        return "de.unijena.cheminf.deglycosylation.tools.DynamicSMILESFileReader";
    }

    @Override
    public String getWriterClassName() {
        return "org.openscience.cdk.io.SMILESWriter";
    }

    @Override
    public int getSupportedDataFeatures() {
        return getRequiredDataFeatures() | DataFeatures.HAS_GRAPH_REPRESENTATION;
    }

    @Override
    public int getRequiredDataFeatures() {
        return DataFeatures.HAS_ATOM_ELEMENT_SYMBOL;
    }

    @Override
    public String getFormatName() {
        return "SMILES";
    }

    @Override
    public String getPreferredNameExtension() {
        return this.getNameExtensions()[0];
    }

    @Override
    public String[] getNameExtensions() {
        return new String[]{".smi", ".txt", ".csv", ".tsv"};
    }

    @Override
    public String getMIMEType() {
        return "chemical/x-daylight-smiles";
    }

    @Override
    public boolean isXMLBased() {
        return false;
    }
}
