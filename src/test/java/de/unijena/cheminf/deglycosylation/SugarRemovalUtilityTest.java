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
 * TODO:
 * - write docs
 * - update the used COCONUT version to the latest publicly available one
 * - print sugar removal util settings
 * - see to dos at failing tests
 */

import com.mongodb.MongoClientSettings;
import com.mongodb.MongoTimeoutException;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoClients;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.MongoDatabase;
import org.bson.Document;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * TODO: Add doc
 */
public class SugarRemovalUtilityTest {
    /**
     * TODO
     * @throws Exception
     */
    @Test
    public void testCircularSugarRemoval() throws Exception {
        //TODO: Use a map?
        String[] tmpSmilesArray = {"CC(N)C(=O)NC(CCC(N)=O)C(=O)NOC1OC(O)C(O)C(O)C1O", //CHEMBL56258
                "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O", //CHEMBL168422
                 "OC1OC(O)C(O)C1OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1OC1C(O)C(O)C(O)OC(O)C1O"}; //own creation
        String[] tmpDeglycosylatedSmilesArray = {"O=C(N)CCC(NC(=O)C(N)C)C(=O)NO",
                "O=C(OC1C(O)C(OC(CO)C1O)C=2C(O)=CC(O)=CC2CO)C=CC=CCC(O)C=CC=CCCCCC",
                "OC1C(OCCCCCCCCCCC)OC(OCCCCCCCCCCCCCCCCC)C(O)C1O"};
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical);
        IAtomContainer tmpOriginalMolecule;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(true);
        tmpSugarRemovalUtil.setStructuresToKeepMode(SugarRemovalUtility.StructuresToKeepMode.HEAVY_ATOM_COUNT);
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(5);
        tmpSugarRemovalUtil.setIncludeNrOfAttachedOxygens(true);
        tmpSugarRemovalUtil.setAttachedOxygensToAtomsInRingRatioThreshold(0.5);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(false);
        tmpSugarRemovalUtil.setPropertyOfSugarContainingMolecules(true);
        for (int i = 0; i < tmpSmilesArray.length; i++) {
            tmpOriginalMolecule = tmpSmiPar.parseSmiles(tmpSmilesArray[i]);
            tmpOriginalMolecule = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, false);
            String tmpSmilesAfterDeglycosylation = tmpSmiGen.create(tmpOriginalMolecule);
            System.out.println(tmpSmilesAfterDeglycosylation);
            Assert.assertEquals(tmpDeglycosylatedSmilesArray[i], tmpSmilesAfterDeglycosylation);
        }
    }

    /**
     * @throws Exception
     */
    @Ignore
    @Test
    public void visualInspectionOfCoconutTest() throws Exception {
        MongoClientSettings.Builder tmpBuilder = MongoClientSettings.builder();
        ServerAddress tmpAddress = new ServerAddress("localhost", 27017);
        tmpBuilder.applyToClusterSettings(builder -> builder.hosts(Collections.singletonList(tmpAddress)));
        MongoClientSettings tmpSettings = tmpBuilder.build();
        MongoClient tmpMongoClient = MongoClients.create(tmpSettings);
        String tmpCollectionName = "COCONUTfebruary20";
        MongoDatabase tmpDatabase = tmpMongoClient.getDatabase(tmpCollectionName);
        String tmpDatabaseName = "uniqueNaturalProduct";
        MongoCollection<Document> tmpCollection = tmpDatabase.getCollection(tmpDatabaseName);
        MongoCursor<Document> tmpCursor = null;
        Logger tmpLogger = Logger.getLogger(SugarRemovalUtilityTest.class.getName());
        try {
            tmpCursor = tmpCollection.find().iterator();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            tmpLogger.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Connection to MongoDB successful.");
        System.out.println("Collection " + tmpCollectionName + " in database " + tmpDatabaseName + " is loaded.");
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        String tmpOutputFolderPath = (new File("SugarRemovalUtilityTest_Output")).getAbsolutePath() + File.separator;
        System.out.println("Output directory: " + tmpOutputFolderPath);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputFolderPath + "Log.txt");
        } catch (IOException anIOException) {
            tmpLogger.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger.getLogger("").addHandler(tmpLogFileHandler);
        Logger.getLogger("").setLevel(Level.ALL);
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(true);
        tmpSugarRemovalUtil.setStructuresToKeepMode(SugarRemovalUtility.StructuresToKeepMode.HEAVY_ATOM_COUNT);
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(5);
        tmpSugarRemovalUtil.setIncludeNrOfAttachedOxygens(true);
        tmpSugarRemovalUtil.setAttachedOxygensToAtomsInRingRatioThreshold(0.5);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(false);
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        int tmpContainsLinearSugarsCounter = 0;
        int tmpContainsCircularSugarsCounter = 0;
        int tmpBasicallyASugarCounter = 0;
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
                tmpID = tmpCurrentDoc.getString("coconut_id");
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                tmpMolecule = tmpSugarRemovalUtil.removeAllSugars(tmpMolecule, true);
                if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    //tmpDepictionGenerator.depict(tmpSmiPar.parseSmiles(tmpSmilesCode)).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                    tmpSugarContainingMoleculesCounter++;
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsCircularSugarsCounter++;
                    }
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsLinearSugarsCounter++;
                    }
                    if (tmpMolecule.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                    }
                }
            } catch (Exception anException) {
                tmpLogger.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpSugarContainingMoleculesCounter);
        System.out.println("Linear-sugar-containing molecules counter: " + tmpContainsLinearSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpContainsCircularSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpCursor.close();
    }

    /**
     *
     *
     * @throws Exception
     */
    @Test
    public void specificTest1() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC2=CC(=CC(OC3OC(CO)C(O)C(O)C3O)=C2C4=C1CCC4)C"); //CNP0000001
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //A simple example, the sugar is not terminal and therefore removed; The resulting disconnected CH3OH is too
        // small to keep and gets cleared away
        Assert.assertEquals("O=C1OC=2C=C(C=C(O)C2C3=C1CCC3)C", tmpSmilesCode);
    }

    /**
     *
     *
     * @throws Exception
     */
    @Test
    public void specificTest2() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4"); //CNP0000012
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The sugar ring is not terminal and should not be removed, so the molecule remains unchanged
        Assert.assertEquals("O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that all sugars are removed, the sugar ring is removed and an unconnected structure remains
        Assert.assertEquals("O=C(O)C=CC1=CC=C(O)C=C1.O=C(OCC1(O)C(O)CC2C(=COC(OC(=O)CC(C)C)C21)CO)C", tmpSmilesCode);
    }

    /**
     *
     *
     * @throws Exception
     */
    @Test
    public void specificTest3() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O"); //CNP0000006
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing is removed, the sugar is terminal because the two phosphate groups are big enough to keep
        Assert.assertEquals("O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O", tmpSmilesCode);
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(6);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, one of the phosphate groups is removed because it has only 5 heavy atoms and therefore, the sugar is
        // no longer terminal and also removed
        Assert.assertEquals("O=P(O)(O)OC", tmpSmilesCode);
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(7);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, both phosphate groups are removed because they are too small and nothing remains of the molecule
        Assert.assertEquals("", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        //back to default
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(5);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, also non-terminal sugars are removed, which leaves two unconnected phosphate groups in this case
        Assert.assertEquals("O=P(O)(O)O.O=P(O)(O)OC", tmpSmilesCode);
        tmpSugarRemovalUtil.setAttachedOxygensToAtomsInRingRatioThreshold(0.7);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, the sugar ring does not have enough oxygen atom attached to be classified as a sugar and be removed
        Assert.assertEquals("O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O", tmpSmilesCode);
    }

    /**
     *
     *
     * @throws Exception
     */
    @Test
    public void specificTest4() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC"); //CNP0000030
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Although there is a match for a linear sugar, it should not be removed
        Assert.assertEquals("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The linear sugar match is not terminal!
        Assert.assertEquals("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //With non-terminal sugars and linear sugars in rings removed, the structure gets disconnected
        Assert.assertEquals("OCCC1CCC2(C=C3C=CCC4C(C=CC(C34)C2C1)CCCC(C)(CC5=CC=C(N)[NH+]=C5)CCC=CC)NCC.C=1C=C2C3=C(NC(NC)CC3CCC2CCC)[NH+]1", tmpSmilesCode);
    }

    /**
     *
     *
     * @throws Exception
     */
    @Test
    public void specificTest5() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //TODO: In this molecule, a CH3 group is removed which should not happen. Cause: There are two overlapping matches of linear sugar patterns
        // and the second match gets reduced to the CH3 group because the overlapping atoms are removed from it. Solution: Combine overlapping linear structures?
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=CCC12OC(OC1=O)C3(C)C(C)CCC3(C)C2"); //CNP0000023
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        System.out.println(tmpSugarRemovalUtil.hasLinearSugars(tmpOriginalMolecule));
        //Nothing should be removed here although there is a match for the linear sugar patterns
        Assert.assertEquals("O=CCC12OC(OC1=O)C3(C)C(C)CCC3(C)C2", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest6() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //TODO: In this molecule, some OH and CH3OH groups are removed which should not happen. Cause: There are two overlapping matches of linear sugar patterns
        // and the second match gets reduced to the small groups because the overlapping atoms are removed from it. Solution: Combine overlapping linear structures?
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O"); //CNP0000082
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is a match for the linear sugar patterns
        Assert.assertEquals("O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest7() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //TODO: In this molecule, a hydroxy group is removed which should not happen. Cause: There are two overlapping matches of linear sugar patterns
        // and the second match gets reduced to the hydroxy group because the overlapping atoms are removed from it. Solution: Combine overlapping linear structures?
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=N[CH+]C5)C(C)C)C(O)C1O"); //CNP0000233
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is a match for the linear sugar patterns
        Assert.assertEquals("O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=N[CH+]C5)C(C)C)C(O)C1O", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest8() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //TODO: In this molecule, a CH3 group is removed which should not happen. Cause: There are two overlapping matches of linear sugar patterns
        // and the second match gets reduced to the CH3 group because the overlapping atoms are removed from it. Solution: Combine overlapping linear structures?
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC(C2=COC=C2)CC3(C)C1CCC4(C)C3C5OC(=O)C4(O)C=C5"); //CNP0000304
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is a match for the linear sugar patterns
        Assert.assertEquals("O=C1OC(C2=COC=C2)CC3(C)C1CCC4(C)C3C5OC(=O)C4(O)C=C5", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest9() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC2CC3(OC4(O)C(CC5(OC45C(=O)OC)CCCCCCCCCCCCCCCC)C2(O3)C1)CCCCCCCCCCCCCCCC"); //CNP0000445
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is a match for the linear sugar patterns
        Assert.assertEquals("O=C1OC2CC3(OC4(O)C(CC5(OC45C(=O)OC)CCCCCCCCCCCCCCCC)C2(O3)C1)CCCCCCCCCCCCCCCC", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest10() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //TODO: In this molecule, a CH3OH group is removed which should not happen. Cause: There are two overlapping matches of linear sugar patterns
        // and the second match gets reduced to the CH3OH group because the overlapping atoms are removed from it. Solution: Combine overlapping linear structures?
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO"); //CNP0000509
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is a match for the linear sugar patterns
        Assert.assertEquals("O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest11() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(OC1OC(CO)C(O)C(O)C1O)(C)CC(=O)OCC=CC2=CC(OC)=C(OC3OC(CO)C(O)C(O)C3O)C(OC)=C2"); //CNP0000920
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The two circular sugar moieties and one connected linear sugar are removed
        Assert.assertEquals("OC1=C(OC)C=C(C=CC)C=C1OC", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest12() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1=CC=C(OC2OC(CO)C(O)C(O)C2O)C=C1"); //CNP0001189
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //One circular and one linear sugar moiety are removed
        Assert.assertEquals("OC1=CC=C(C=C1)C", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest13() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(COC(=O)C(CC(=O)O)C(OCC3C(=C)CCC4C(C)(C)CCCC34C)C(=O)OC)CCCC12C)C(=O)OC"); //CNP0001863
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //One linear sugar moiety is removed, the other one is not because it is terminal; no other changes are done to the molecule
        Assert.assertEquals("O=C(O)CC(C(=O)OCC1(C)CCCC2(C)C(C(=C)CCC12)C)C(OCC3C(=C)CCC4C(C)(C)CCCC34C)C(=O)OC", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest14() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2OC%10OC(CO)C(O)C(O)C%10O)CO)C(O)C1O"); //CNP0002871
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("OC1CC2CCC3C(CCC4(C)C3CC5OC6(OCC(C)CC6)C(C)C54)C2(C)CC1O", tmpSmilesCode);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC4C5=CCC6C7(C)CCC(O)C(C(=O)OC8OC(CO)C(O)C(O)C8O)(C)C7CCC6(C)C5(C)CC3)C(O)C(OC9OC(CO)C(O)C(O)C9O)C2O)C(OC%10OC(CO)C(O)C(O)C%10O)C(O)C1O"); //CNP0005247
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(O)C1(C)C(O)CCC2(C)C1CCC3(C)C2CC=C4C5CC(C)(C)CCC5(C(=O)O)CCC43C", tmpSmilesCode);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OC1C(O)C(OC2C3=C(O)C(=CC=C3OC2C(=C)CO)C(=O)C)OC(CO)C1O"); //CNP0032326
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(C1=CC=C2OC(C(=C)CO)C(O)C2=C1O)C", tmpSmilesCode);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"); //CNP0031401
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C1C=C(OC=2C=C(O)C=C(O)C12)C=3C=CC(O)=C(O)C3", tmpSmilesCode);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2OC3OCC(O)C(OC4OCC(O)C(O)C4O)C3O)OC5C(OC(=O)C67CCC(C)(C)CC7C8=CCC9C%10(C)CC(O)C(OC%11OC(CO)C(O)C(O)C%11O)C(C(=O)O)(C)C%10CCC9(C)C8(CO)CC6)OC(C)C(OC(=O)C=CC%12=CC(OC)=C(OC)C(OC)=C%12)C5OC%13OC(C)C(O)C(O)C%13O)C1O"); //CNP0028122
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(OC1C(O)C(O)C(OC(=O)C23CCC(C)(C)CC3C4=CCC5C6(C)CC(O)C(O)C(C(=O)O)(C)C6CCC5(C)C4(CO)CC2)OC1C)C=CC7=CC(OC)=C(OC)C(OC)=C7", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest15() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C1OC(O)C(O)C(O)C1O"); //CNP0002919
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //This molecule is a sugar and should be completely removed
        Assert.assertEquals("", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest16() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(C)CCCC12C)C(=O)O"); //CNP0003329
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The linear sugar moieties are removed
        Assert.assertEquals("C=C1CCC2C(C)(C)CCCC2(C)C1C", tmpSmilesCode);
        //another example CNP0031156
    }

    /**
     *
     */
    @Test
    public void specificTest17() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C=C(OC=2C1=CC3=C(OC(C)(C)C(OOCC(O)C(O)C(O)C(O)CO)C3)C2[N+]=4C=C5N=CC=C5C4CC)C"); //CNP0007654
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The linear sugar moieties are removed
        Assert.assertEquals("O=C1C=C(OC=2C1=CC3=C(OC(C)(C)C(O)C3)C2[N+]=4C=C5N=CC=C5C4CC)C", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest18() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4"); //CNP0032817
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The sugar moiety is NOT connected to the core structure via a glycosidic bond, so it is not removed
        Assert.assertEquals("O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that this setting is changed, the sugar moiety is removed
        Assert.assertEquals("O=COC=1C=C2OC(=CC(=O)C2=C(O)C1O)C=3C=CC(O)=CC3", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        //another examples for the same thing:
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"); //CNP0031401
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The circular sugar moiety is NOT connected to the core structure via a glycosidic bond, so it is not removed
        Assert.assertEquals("O=C1C=C(OC=2C1=C(O)C=C(O)C2C3OC(CO)C(O)C(O)C3O)C=4C=CC(O)=C(O)C4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that this setting is changed, the sugar moiety is removed
        Assert.assertEquals("O=C1C=C(OC=2C=C(O)C=C(O)C12)C=3C=CC(O)=C(O)C3", tmpSmilesCode);
    }

    /**
     * This molecule caused an exception before because a linear sugar candidate got removed while too small structures
     * were cleared away.
     */
    @Test
    public void specificTest19() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC"); //CNP0000913
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeAllSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Only the aliphatic chain remains
        Assert.assertEquals("CCCCCCCCCCC", tmpSmilesCode);
    }

    /**
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void mappingsExperiment() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer target = tmpSmiPar.parseSmiles("CCCCC");
        IAtomContainer query = tmpSmiPar.parseSmiles("CC");
        DfPattern tmpPattern = DfPattern.findSubstructure(query);
        Mappings tmpMappings = tmpPattern.matchAll(target);
        for (IAtomContainer tmpMap : tmpMappings.toSubstructures()) {
            System.out.println(tmpSmiGen.create(tmpMap));
        }
        System.out.println(tmpMappings.countUnique());
    }
}
