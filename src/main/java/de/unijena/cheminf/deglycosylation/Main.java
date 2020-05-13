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

/**
 * TODO
 */
public class Main {
    /**
     * TODO
     * See Readers in COCONUT
     * Output: Meditate about it (the detailed structure of the output file) CSV?
     * Should include the SMILES of the original structure and the removed moieties
     *
     * arg 0: SMILES code or path to SMILES file, SDF, or molfile [String]
     * arg 1: param to indicate whether circular, linear, or both types of sugars should be removed [int?]
     * the other params should be optional (but if one is present, all should be present?)
     * arg 2: detect glycosidic bonds [boolean]
     * arg 3: remove only terminal [boolean]
     *      -> this requires checking all the molecules for being connected before processing them!
     * arg 4: structure to keep mode [int?]
     * arg 5: structure to keep mode threshold [int]
     * arg 6: include number of attached oxygens [boolean]
     * arg 7: ratio attached oxygens to atoms in ring threshold [double]
     * arg 8: remove linear sugars in ring [boolean]
     * arg 9: linear sugar candidate min size [int]
     * arg 10: linear sugar candidate max size [int]
     */
    public static void main(String[] args) {

    }
}
