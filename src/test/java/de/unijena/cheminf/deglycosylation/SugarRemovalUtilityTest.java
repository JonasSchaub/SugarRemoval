/**
 * TODO: Add license header (switch to MIT)
 */
package de.unijena.cheminf.deglycosylation;

/**
 * TODO:
 * - Add test for removal of linear sugars!
 * - Add more examples
 * - play around with settings
 * - visually inspect results of coconut test, also to get more test cases
 * - write docs
 */

import com.mongodb.MongoClientSettings;
import com.mongodb.MongoTimeoutException;
import com.mongodb.ServerAddress;
import com.mongodb.client.*;
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
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Objects;
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
        //IAtomContainer tmpDeglycosylatedMolecule;
        //UniversalIsomorphismTester tmpUnivIsomorphTester = new UniversalIsomorphismTester();
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
            //tmpDeglycosylatedMolecule = tmpSmiPar.parseSmiles(tmpDeglycosylatedSmilesArray[i]);
            tmpOriginalMolecule = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, false);
            String tmpSmilesAfterDeglycosylation = tmpSmiGen.create(tmpOriginalMolecule);
            System.out.println(tmpSmilesAfterDeglycosylation);
            Assert.assertEquals(tmpDeglycosylatedSmilesArray[i], tmpSmilesAfterDeglycosylation);
            //Note: for an unknown reason, this does not work all the time...
            //Assert.assertTrue(tmpUnivIsomorphTester.isIsomorph(tmpOriginalMolecule, tmpDeglycosylatedMolecule));
        }
    }

    /**
     *
     * @throws Exception
     */
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
