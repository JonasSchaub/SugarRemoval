/**
 * TODO: Add license header (switch to MIT?)
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
                tmpMolecule = tmpSugarRemovalUtil.removeAllSugars(tmpMolecule, false);
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
