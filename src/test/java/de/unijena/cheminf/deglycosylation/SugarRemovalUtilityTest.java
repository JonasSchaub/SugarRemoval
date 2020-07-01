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
 * - linear sugar detection: make some inspections about the non-default settings
 * - include tests for static methods
 * - test the protected routines
 * - remove some of the 'experiments' in the end
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
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.*;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.stream.Collectors;

/**
 * TODO: Add doc
 */
public class SugarRemovalUtilityTest extends SugarRemovalUtility {
    //<editor-fold desc="Tests involving databases">
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
        String tmpCollectionName = "COCONUTmay";
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
        String tmpOutputFolderPath = (new File("SugarRemovalUtilityTest_Output")).getAbsolutePath() + File.separator
                + "coconut_inspection_test" + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
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
        Logger.getLogger("").setLevel(Level.WARNING);
        //Done for reproducibility
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        System.out.println(tmpSugarRemovalUtil.areOnlyCircularSugarsWithOGlycosidicBondDetected());
        System.out.println(tmpSugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeSetting());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeThresholdSetting());
        System.out.println(tmpSugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected());
        System.out.println(tmpSugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        System.out.println(tmpSugarRemovalUtil.areLinearSugarsInRingsDetected());
        System.out.println(tmpSugarRemovalUtil.arePropertiesAddedToSugarContainingMolecules());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMinSizeSetting());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting());
        System.out.println(tmpSugarRemovalUtil.areLinearAcidicSugarsDetected());
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
                IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    //tmpDepictionGenerator.depict(tmpDeglycosylatedClone).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                    tmpSugarContainingMoleculesCounter++;
                    if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsCircularSugarsCounter++;
                    }
                    if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                        //tmpDepictionGenerator.depict(tmpDeglycosylatedClone).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                        tmpContainsLinearSugarsCounter++;
                    }
                    /*if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true
                    && (boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                        //tmpDepictionGenerator.depict(tmpDeglycosylatedClone).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                    }*/
                    if (tmpDeglycosylatedClone.isEmpty()) {
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
     * TODO
     */
    @Ignore
    @Test
    public void CoconutStatsTest() throws Exception {

        //<editor-fold desc="Connection to MongoDB">
        MongoClientSettings.Builder tmpBuilder = MongoClientSettings.builder();
        ServerAddress tmpAddress = new ServerAddress("localhost", 27017);
        tmpBuilder.applyToClusterSettings(builder -> builder.hosts(Collections.singletonList(tmpAddress)));
        MongoClientSettings tmpSettings = tmpBuilder.build();
        MongoClient tmpMongoClient = MongoClients.create(tmpSettings);
        String tmpCollectionName = "COCONUTmay";
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
        String tmpOutputFolderPath = (new File("SugarRemovalUtilityTest_Output")).getAbsolutePath() + File.separator
                + "coconut_stats_test" + File.separator;
        //</editor-fold>

        //<editor-fold desc="Output folder and logging setup">
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
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
        Logger.getLogger("").setLevel(Level.WARNING);
        //</editor-fold>

        NumberFormat tmpRatioOutputFormat = NumberFormat.getInstance(Locale.US);
        tmpRatioOutputFormat.setMaximumFractionDigits(1);
        tmpRatioOutputFormat.setRoundingMode(RoundingMode.DOWN);

        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");

        //Done for reproducibility
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();

        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;

        //<editor-fold desc="Counter and map definitions">
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        List<String> tmpHasAnyTypeOfSugarsCNPs = new ArrayList<>(50000);
        int tmpHasNoSugarsCounter = 0;
        //List<String> tmpHasNoSugarsCNPs = new ArrayList<>(400000);

        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasNonTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasOnlyTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasOnlyNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasTerminalAndNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalCircularSugarsCNPs = new ArrayList(5000);
        int tmpHasOnlyCircularSugarsCounter = 0;
        List<String> tmpHasOnlyCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasGlycosidicBondCounter = 0;
        List<String> tmpHasGlycosidicBondCNPs = new ArrayList(50000);
        int tmpHasGlycosidicBondOnTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnTerminalSugarCNPs = new ArrayList(50000);
        int tmpHasGlycosidicBondOnNonTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnNonTerminalSugarCNPs = new ArrayList(50000);
        int tmpGlycosidicBondExemptionCounter = 0;
        List<String> tmpGlycosidicBondExemptionCNPs = new ArrayList(120);

        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasNonTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasOnlyTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasOnlyNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasTerminalAndNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalLinearSugarsCNPs = new ArrayList(10);
        int tmpHasOnlyLinearSugarsCounter = 0;
        List<String> tmpHasOnlyLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasLinearSugarsInRingCounter = 0;
        List<String> tmpHasLinearSugarsInRingCNPs = new ArrayList(500);

        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList(2000);

        int tmpBasicallyASugarCounter = 0; //number of molecules that are basically sugars, circular or linear, polymer or single unit
        List<String> tmpBasicallyASugarCNPs = new ArrayList(2000);
        int tmpBasicallyACircularSugarCounter = 0;
        List<String> tmpBasicallyACircularSugarCNPs = new ArrayList(2000);
        int tmpBasicallyALinearSugarCounter = 0;
        List<String> tmpBasicallyALinearSugarCNPs = new ArrayList(2000);
        int tmpBasicallyASingleSugarUnitCounter = 0; //circular or linear
        List<String> tmpBasicallyASingleSugarUnitCNPs = new ArrayList(1000);
        int tmpBasicallyASingleCircularSugarCounter = 0;
        List<String> tmpBasicallyASingleCircularSugarCNPs = new ArrayList(1000);
        int tmpBasicallyASingleLinearSugarCounter = 0;
        List<String> tmpBasicallyASingleLinearSugarCNPs = new ArrayList(1000);

        int tmpCircularSugarMoietiesCounter = 0;
        int tmpTerminalCircularSugarMoietiesCounter = 0;
        int tmpNonTerminalCircularSugarMoietiesCounter = 0;
        int tmpCircularSugarMoietiesWithGlycosidicBondCounter = 0; //the rest have no glycosidic bond
        int tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter = 0; //the rest are non-terminal

        int tmpLinearSugarMoietiesCounter = 0;
        int tmpTerminalLinearSugarMoietiesCounter = 0;
        int tmpNonTerminalLinearSugarMoietiesCounter = 0;
        int tmpLinearSugarMoietiesInRingsCounter = 0;

        HashMap<Integer, Integer> tmpFrequenciesOfSizesOfCircularSugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap = new HashMap(10, 0.9f);

        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManySugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap = new HashMap(10, 0.9f);

        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap = new HashMap(10, 0.9f);
        HashMap<String, Integer> tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap = new HashMap(10, 0.9f);
        int tmpUnexpectedRingSizeCounter = 0;

        int tmpLinearSugarsDetectedInCircularSugarsCounter = 0;
        //</editor-fold>

        //<editor-fold desc="Iteration of molecules">
        while (tmpCursor.hasNext()) {
            try {

                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString("inchikey"); //coconut_id or inchikey
                tmpSmilesCode = tmpCurrentDoc.getString("smiles"); //smiles or clean_smiles
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);

                //using default settings where nothing else is specified
                boolean tmpHasAnyTypeOfSugar = tmpSugarRemovalUtil.hasCircularAndOrLinearSugars(tmpMolecule);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);

                if (tmpHasAnyTypeOfSugar) {

                    tmpHasAnyTypeOfSugarsCounter++;
                    tmpHasAnyTypeOfSugarsCNPs.add(tmpID);

                    int tmpNumberOfCircularAndLinearSugarMoieties = tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpMolecule);
                    if (!tmpHowManyMoleculesHaveHowManySugarMoietiesMap.containsKey(tmpNumberOfCircularAndLinearSugarMoieties)) {
                        tmpHowManyMoleculesHaveHowManySugarMoietiesMap.put(tmpNumberOfCircularAndLinearSugarMoieties, 1);
                    } else {
                        Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(tmpNumberOfCircularAndLinearSugarMoieties);
                        tmpHowManyMoleculesHaveHowManySugarMoietiesMap.put(tmpNumberOfCircularAndLinearSugarMoieties, tmpCurrentListValue + 1);
                    }

                    if (tmpHasAnyCircularSugar) {

                        //<editor-fold desc="Analysis of molecules having circular sugars">
                        tmpHasCircularSugarsCounter++;
                        tmpHasCircularSugarsCNPs.add(tmpID);

                        if (!tmpHasAnyLinearSugar) {
                            tmpHasOnlyCircularSugarsCounter++;
                            tmpHasOnlyCircularSugarsCNPs.add(tmpID);
                        }

                        //terminal and non-terminal, having a glycosidic bond or not (see default settings)
                        List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());

                        int tmpNumberOfCircularSugarMoieties;
                        int tmpNumberOfTerminalCircularSugarMoieties;
                        int tmpNumberOfNonTerminalCircularSugarMoieties;
                        int tmpNumberOfGlycosidicBonds;
                        int tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                        int tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond;

                        tmpNumberOfCircularSugarMoieties = tmpCircularSugarCandidatesList.size();

                        if (!tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.containsKey(tmpNumberOfCircularSugarMoieties)) {
                            tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.put(tmpNumberOfCircularSugarMoieties, 1);
                        } else {
                            Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(tmpNumberOfCircularSugarMoieties);
                            tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.put(tmpNumberOfCircularSugarMoieties, tmpCurrentValue + 1);
                        }

                        //the first is the overall counter, the second one is specific for this molecule
                        tmpCircularSugarMoietiesCounter += tmpNumberOfCircularSugarMoieties;

                        for (IAtomContainer tmpCircularSugarCandidate : tmpCircularSugarCandidatesList) {
                            int tmpCandidateSize = AtomContainerManipulator.getHeavyAtoms(tmpCircularSugarCandidate).size();
                            if (!tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.containsKey(tmpCandidateSize)) {
                                tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.put(tmpCandidateSize, 1);
                            } else {
                                Integer tmpCurrentCount = tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(tmpCandidateSize);
                                tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.put(tmpCandidateSize, tmpCurrentCount + 1);
                            }
                        }

                        //note: circular moieties that become terminal after removal of a linear moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalCircularSugarMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        //-1 for the deglycosylated core at the beginning of the list
                        tmpNumberOfTerminalCircularSugarMoieties = tmpRemovedTerminalCircularSugarMoieties.size() - 1 ;
                        tmpNumberOfNonTerminalCircularSugarMoieties = tmpNumberOfCircularSugarMoieties - tmpNumberOfTerminalCircularSugarMoieties;
                        Assert.assertTrue(tmpNumberOfNonTerminalCircularSugarMoieties >= 0);

                        //leaving default! Now, only circular sugars having glycosidic bonds are in the candidates and removed moieties
                        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
                        boolean tmpMoleculeQualifiesForExemption = tmpSugarRemovalUtil.isQualifiedForGlycosidicBondExemption(tmpMolecule.clone());
                        if (tmpMoleculeQualifiesForExemption) {
                            //note: these molecules are basically made up of one circular sugar without a glycosidic bond,
                            // so that the whole molecule is empty after the removal of this sugar moiety. This is an
                            // exemption implemented in the compilation of circular sugar moieties with detection of glycosidic bonds.
                            tmpGlycosidicBondExemptionCounter++;
                            tmpGlycosidicBondExemptionCNPs.add(tmpID);
                            tmpNumberOfGlycosidicBonds = 0;
                            tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond = 0;
                            tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond = 0;
                        } else {
                            List<IAtomContainer> tmpCircularSugarCandidatesWithGlycosidicBondsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                            tmpNumberOfGlycosidicBonds = tmpCircularSugarCandidatesWithGlycosidicBondsList.size();
                            List<IAtomContainer> tmpRemovedTerminalCircularMoietiesWithGlycosidicBond = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                            tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond = tmpRemovedTerminalCircularMoietiesWithGlycosidicBond.size() - 1;
                            tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond = tmpNumberOfGlycosidicBonds - tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                            Assert.assertTrue(tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond >= 0);
                        }
                        //back to default!
                        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);

                        if (tmpNumberOfTerminalCircularSugarMoieties > 0) {

                            tmpHasTerminalCircularSugarsCounter++;
                            tmpHasTerminalCircularSugarsCNPs.add(tmpID);
                            tmpTerminalCircularSugarMoietiesCounter += tmpNumberOfTerminalCircularSugarMoieties;

                            if (tmpNumberOfNonTerminalCircularSugarMoieties == 0) {
                                tmpHasOnlyTerminalCircularSugarsCounter++;
                                tmpHasOnlyTerminalCircularSugarsCNPs.add(tmpID);
                            }

                        }

                        if (tmpNumberOfNonTerminalCircularSugarMoieties > 0) {

                            tmpHasNonTerminalCircularSugarsCounter++;
                            tmpHasNonTerminalCircularSugarsCNPs.add(tmpID);
                            tmpNonTerminalCircularSugarMoietiesCounter += tmpNumberOfNonTerminalCircularSugarMoieties;

                            if (tmpNumberOfTerminalCircularSugarMoieties == 0) {
                                tmpHasOnlyNonTerminalCircularSugarsCounter++;
                                tmpHasOnlyNonTerminalCircularSugarsCNPs.add(tmpID);
                            }

                        }

                        if (tmpNumberOfTerminalCircularSugarMoieties > 0 && tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                            tmpHasTerminalAndNonTerminalCircularSugarsCounter++;
                            tmpHasTerminalAndNonTerminalCircularSugarsCNPs.add(tmpID);
                        }

                        if (tmpNumberOfGlycosidicBonds > 0) {

                            tmpHasGlycosidicBondCounter++;
                            tmpHasGlycosidicBondCNPs.add(tmpID);

                            //the first is the overall counter, the second one is specific for this molecule
                            tmpCircularSugarMoietiesWithGlycosidicBondCounter += tmpNumberOfGlycosidicBonds;

                            if (!tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.containsKey(tmpNumberOfGlycosidicBonds)) {
                                tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.put(tmpNumberOfGlycosidicBonds, 1);
                            } else {
                                Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(tmpNumberOfGlycosidicBonds);
                                tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.put(tmpNumberOfGlycosidicBonds, tmpCurrentListValue + 1);
                            }

                            if (tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond > 0) {
                                tmpHasGlycosidicBondOnTerminalSugarCounter++;
                                tmpHasGlycosidicBondOnTerminalSugarCNPs.add(tmpID);
                                //the first is the overall counter, the second one is specific for this molecule
                                tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter += tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                            }

                            if (tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond > 0) {
                                tmpHasGlycosidicBondOnNonTerminalSugarCounter++;
                                tmpHasGlycosidicBondOnNonTerminalSugarCNPs.add(tmpID);
                            }

                        }
                        //</editor-fold>

                    }

                    //<editor-fold desc="Analysis of molecules having linear sugars">
                    if (tmpHasAnyLinearSugar) {

                        tmpHasLinearSugarsCounter++;
                        tmpHasLinearSugarsCNPs.add(tmpID);

                        if (!tmpHasAnyCircularSugar) {
                            tmpHasOnlyLinearSugarsCounter++;
                            tmpHasOnlyLinearSugarsCNPs.add(tmpID);
                        }

                        //terminal and non-terminal
                        List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule.clone());

                        int tmpNumberOfLinearSugarMoieties;
                        int tmpNumberOfTerminalLinearSugarMoieties;
                        int tmpNumberOfNonTerminalLinearSugarMoieties;

                        tmpNumberOfLinearSugarMoieties = tmpLinearSugarCandidatesList.size();

                        //the first is the overall counter, the second one is specific for this molecule
                        tmpLinearSugarMoietiesCounter += tmpNumberOfLinearSugarMoieties;

                        if (!tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.containsKey(tmpNumberOfLinearSugarMoieties)) {
                            tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.put(tmpNumberOfLinearSugarMoieties, 1);
                        } else {
                            Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(tmpNumberOfLinearSugarMoieties);
                            tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.put(tmpNumberOfLinearSugarMoieties, tmpCurrentValue + 1);
                        }

                        for (IAtomContainer tmpLinearSugarCandidate : tmpLinearSugarCandidatesList) {

                            int tmpCandidateSize = AtomContainerManipulator.getHeavyAtoms(tmpLinearSugarCandidate).size();
                            if (!tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.containsKey(tmpCandidateSize)) {
                                tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.put(tmpCandidateSize, 1);
                            } else {
                                Integer tmpCurrentCount = tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(tmpCandidateSize);
                                tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.put(tmpCandidateSize, tmpCurrentCount + 1);
                            }

                            int tmpCarbonCount = 0;
                            for (IAtom tmpAtom : tmpLinearSugarCandidate.atoms()) {
                                String tmpSymbol = tmpAtom.getSymbol();
                                if (tmpSymbol.equals("C")) {
                                    tmpCarbonCount++;
                                }
                            }
                            if (!tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.containsKey(tmpCarbonCount)) {
                                tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.put(tmpCarbonCount, 1);
                            } else {
                                Integer tmpCurrentCount = tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(tmpCarbonCount);
                                tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.put(tmpCarbonCount, tmpCurrentCount + 1);
                            }

                        }

                        //note:linear moieties that become terminal after removal of a circular moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalLinearMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                        //-1 for the deglycosylated core at the beginning of the list
                        tmpNumberOfTerminalLinearSugarMoieties = tmpRemovedTerminalLinearMoieties.size() - 1 ;
                        tmpNumberOfNonTerminalLinearSugarMoieties = tmpNumberOfLinearSugarMoieties - tmpNumberOfTerminalLinearSugarMoieties;
                        Assert.assertTrue(tmpNumberOfNonTerminalLinearSugarMoieties >= 0);

                        if (tmpNumberOfTerminalLinearSugarMoieties > 0) {

                            tmpHasTerminalLinearSugarsCounter++;
                            tmpHasTerminalLinearSugarsCNPs.add(tmpID);
                            tmpTerminalLinearSugarMoietiesCounter += tmpNumberOfTerminalLinearSugarMoieties;

                            if (tmpNumberOfNonTerminalLinearSugarMoieties == 0) {
                                tmpHasOnlyTerminalLinearSugarsCounter++;
                                tmpHasOnlyTerminalLinearSugarsCNPs.add(tmpID);
                            }

                        }

                        if (tmpNumberOfNonTerminalLinearSugarMoieties > 0) {

                            tmpHasNonTerminalLinearSugarsCounter++;
                            tmpHasNonTerminalLinearSugarsCNPs.add(tmpID);
                            tmpNonTerminalLinearSugarMoietiesCounter += tmpNumberOfNonTerminalLinearSugarMoieties;

                            if (tmpNumberOfTerminalLinearSugarMoieties == 0) {
                                tmpHasOnlyNonTerminalLinearSugarsCounter++;
                                tmpHasOnlyNonTerminalLinearSugarsCNPs.add(tmpID);
                            }

                        }

                        if (tmpNumberOfTerminalLinearSugarMoieties > 0 && tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                            tmpHasTerminalAndNonTerminalLinearSugarsCounter++;
                            tmpHasTerminalAndNonTerminalLinearSugarsCNPs.add(tmpID);
                        }

                    }
                    //</editor-fold>

                    if (tmpHasAnyCircularSugar && tmpHasAnyLinearSugar) {
                        tmpHasCircularAndLinearSugarsCounter++;
                        tmpHasCircularAndLinearSugarsCNPs.add(tmpID);
                    }

                    //<editor-fold desc="Analysis of molecules that are basically sugars">
                    //removes only terminal moieties but that is correct here
                    List<IAtomContainer> tmpDeglycosylatedCloneAndRemovedSugarMoietiesList =
                            tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpMolecule, true);
                    IAtomContainer tmpDeglycosylatedClone = tmpDeglycosylatedCloneAndRemovedSugarMoietiesList.get(0);

                    if (tmpDeglycosylatedClone.isEmpty()) {

                        tmpBasicallyASugarCounter++;
                        tmpBasicallyASugarCNPs.add(tmpID);

                        //note: it is important to count the actually removed moieties here, not the detected ones!
                        // Because there are multiple rounds of detection in the removal if only terminal moieties are removed
                        int tmpNumberOfMoieties = tmpDeglycosylatedCloneAndRemovedSugarMoietiesList.size() - 1;

                        if (tmpNumberOfMoieties == 1) {
                            tmpBasicallyASingleSugarUnitCounter++;
                            tmpBasicallyASingleSugarUnitCNPs.add(tmpID);
                        }

                    }

                    IAtomContainer tmpCircularDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularSugars(tmpMolecule, true);

                    if (tmpCircularDeglycosylatedClone.isEmpty()) {

                        tmpBasicallyACircularSugarCounter++;
                        tmpBasicallyACircularSugarCNPs.add(tmpID);

                        //note: here, it is ok to only count the detected moieties because there is only one round of detection in the removal
                        int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpMolecule.clone());

                        if (tmpNumberOfMoieties == 1) {
                            tmpBasicallyASingleCircularSugarCounter++;
                            tmpBasicallyASingleCircularSugarCNPs.add(tmpID);
                        }

                    }

                    IAtomContainer tmpLinearDeglycosylatedClone = tmpSugarRemovalUtil.removeLinearSugars(tmpMolecule, true);

                    if (tmpLinearDeglycosylatedClone.isEmpty()) {

                        tmpBasicallyALinearSugarCounter++;
                        tmpBasicallyALinearSugarCNPs.add(tmpID);

                        //note: here, it is ok to only count the detected moieties because there is only one round of detection in the removal
                        int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfLinearSugars(tmpMolecule);

                        if (tmpNumberOfMoieties == 1) {
                            tmpBasicallyASingleLinearSugarCounter++;
                            tmpBasicallyASingleLinearSugarCNPs.add(tmpID);
                        }

                    }
                    //</editor-fold>

                } else {
                    tmpHasNoSugarsCounter++;
                    //tmpHasNoSugarsCNPs.add(tmpID);
                }

                //<editor-fold desc="Analysis of all possible sugar rings and their exocyclic oxygen counts">
                //leaving default settings!
                tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(false);
                IAtomContainer tmpMoleculeClone = tmpMolecule.clone();
                List<IAtomContainer> tmpCircularSugarsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMoleculeClone);
                if (tmpCircularSugarsList.size() > 0) {
                    for (IAtomContainer tmpCandidate : tmpCircularSugarsList) {
                        int tmpRingSize = tmpCandidate.getAtomCount();
                        int tmpExocyclicOxygenAtomsCount = this.getExocyclicOxygenAtomCount(tmpCandidate, tmpMoleculeClone);
                        double tmpAttachedOxygensToAtomsInRingRatio =
                                ((double) tmpExocyclicOxygenAtomsCount / (double) tmpRingSize);
                        //note: the ratios are not rounded, the remaining decimals are neglected, which is correct here
                        // because the respective setting is a threshold
                        String tmpRoundedRatio = tmpRatioOutputFormat.format(tmpAttachedOxygensToAtomsInRingRatio);
                        if (!tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.containsKey(tmpRoundedRatio)) {
                            tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.put(tmpRoundedRatio, 1);
                        } else {
                            Integer tmpCurrentCount = tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(tmpRoundedRatio);
                            tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.put(tmpRoundedRatio, tmpCurrentCount + 1);
                        }
                        switch (tmpRingSize) {
                            case 5:
                                if (!tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.containsKey(tmpExocyclicOxygenAtomsCount)) {
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, 1);
                                } else {
                                    Integer tmpCurrentCount = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(tmpExocyclicOxygenAtomsCount);
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, tmpCurrentCount + 1);
                                }
                                break;
                            case 6:
                                if (!tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.containsKey(tmpExocyclicOxygenAtomsCount)) {
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, 1);
                                } else {
                                    Integer tmpCurrentCount = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(tmpExocyclicOxygenAtomsCount);
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, tmpCurrentCount + 1);
                                }
                                break;
                            case 7:
                                if (!tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.containsKey(tmpExocyclicOxygenAtomsCount)) {
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, 1);
                                } else {
                                    Integer tmpCurrentCount = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(tmpExocyclicOxygenAtomsCount);
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, tmpCurrentCount + 1);
                                }
                                break;
                            default:
                                tmpUnexpectedRingSizeCounter++;
                                break;
                        }
                    }
                }
                //back to default settings!
                tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(true);
                //</editor-fold>

                //<editor-fold desc="Analysis of linear sugars in rings">
                //leaving default settings!
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                IAtomContainer tmpNewClone = tmpMolecule.clone();
                List<IAtomContainer> tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                int tmpListSizeWithCandidatesInCycles = tmpLinearCandidates.size();
                if (tmpListSizeWithCandidatesInCycles > 0) {
                    //this.removeSugarCandidatesWithCyclicAtoms(tmpLinearCandidates, tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    int tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidates.size();
                    int tmpNumberOfLinearSugarsInCycles = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    if (tmpNumberOfLinearSugarsInCycles > 0) {
                        tmpHasLinearSugarsInRingCounter++;
                        tmpHasLinearSugarsInRingCNPs.add(tmpID);
                        tmpLinearSugarMoietiesInRingsCounter += tmpNumberOfLinearSugarsInCycles;
                    }
                    //leaving default further!
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
                    tmpSugarRemovalUtil.removeCircularSugars(tmpNewClone, false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpListSizeWithCandidatesInCycles = tmpLinearCandidates.size();
                    //this.removeSugarCandidatesWithCyclicAtoms(tmpLinearCandidates, tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidates.size();
                    int tmpNumberOfLinearSugarsInCyclesWithoutCircularSugars = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    int tmpLinearSugarsThatWereDetectedInCircularSugars = tmpNumberOfLinearSugarsInCycles - tmpNumberOfLinearSugarsInCyclesWithoutCircularSugars;
                    tmpLinearSugarsDetectedInCircularSugarsCounter += tmpLinearSugarsThatWereDetectedInCircularSugars;
                    //back to this default
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
                }
                //back to default settings!
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                //</editor-fold>

            } catch (Exception anException) {
                tmpLogger.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }

        } //end of while()
        //</editor-fold>

        //<editor-fold desc="Printout">
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println();
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println();
        System.out.println("Sugar containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        double tmpPercentage = ((double) tmpHasAnyTypeOfSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain sugars.");
        System.out.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        System.out.println();
        System.out.println("Circular sugar containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpPercentage = ((double) tmpHasCircularSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain circular sugars.");
        System.out.println("Only circular sugar containing molecules counter: " + tmpHasOnlyCircularSugarsCounter);
        System.out.println("The rest contain linear sugars only or both types of sugars (see below)");
        System.out.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalCircularSugarsCounter);
        System.out.println("Only terminal circular sugar containing molecules counter: " + tmpHasOnlyTerminalCircularSugarsCounter);
        System.out.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalCircularSugarsCounter);
        System.out.println("Only non-terminal circular sugar containing molecules counter: " + tmpHasOnlyNonTerminalCircularSugarsCounter);
        System.out.println("Terminal and non-terminal circular sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalCircularSugarsCounter);
        System.out.println("Circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondCounter);
        System.out.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasGlycosidicBondCounter) + " molecules " +
                "only have circular sugar moieties that are NOT attached via a glycosidic bond.");
        System.out.println("Terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnTerminalSugarCounter);
        System.out.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on non-terminal circular sugar moieties.");
        System.out.println("Non-terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnNonTerminalSugarCounter);
        System.out.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnNonTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on terminal circular sugar moieties.");
        int tmpHasBoth = tmpHasGlycosidicBondOnNonTerminalSugarCounter + tmpHasGlycosidicBondOnTerminalSugarCounter - tmpHasGlycosidicBondCounter;
        System.out.println(tmpHasBoth + " molecules have both, terminal and non-terminal sugar moieties attached via a glycosidic bond.");
        System.out.println("Molecules that qualify for the glycosidic bond exemption counter: " + tmpGlycosidicBondExemptionCounter);
        System.out.println();
        System.out.println("Detected circular sugar moieties counter: " + tmpCircularSugarMoietiesCounter);
        System.out.println("Detected terminal circular sugar moieties counter: " + tmpTerminalCircularSugarMoietiesCounter);
        System.out.println("Detected non-terminal circular sugar moieties counter: " + tmpNonTerminalCircularSugarMoietiesCounter);
        System.out.println("Detected circular sugar moieties that have a glycosidic bond counter: " + tmpCircularSugarMoietiesWithGlycosidicBondCounter);
        System.out.println((tmpCircularSugarMoietiesCounter - tmpCircularSugarMoietiesWithGlycosidicBondCounter) + " circular sugar moieties do not have a glycosidic bond.");
        System.out.println("Detected circular sugar moieties that have a glycosidic bond and are terminal counter: " + tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter);
        System.out.println((tmpCircularSugarMoietiesWithGlycosidicBondCounter - tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter)
                + " circular sugar moieties that have a glycosidic bond are non-terminal.");
        System.out.println();
        System.out.println("Linear sugar containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpPercentage = ((double) tmpHasLinearSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain linear sugars.");
        System.out.println("Only linear sugar containing molecules counter: " + tmpHasOnlyLinearSugarsCounter);
        System.out.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        System.out.println("Only terminal linear sugar containing molecules counter: " + tmpHasOnlyTerminalLinearSugarsCounter);
        System.out.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        System.out.println("Only non-terminal linear sugar containing molecules counter: " + tmpHasOnlyNonTerminalLinearSugarsCounter);
        System.out.println("Terminal and non-terminal linear sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalLinearSugarsCounter);
        System.out.println("Linear sugar moieties in rings containing molecules counter: " + tmpHasLinearSugarsInRingCounter);
        System.out.println();
        System.out.println("Detected linear sugar moieties counter: " + tmpLinearSugarMoietiesCounter);
        System.out.println("Detected terminal linear sugar moieties counter: " + tmpTerminalLinearSugarMoietiesCounter);
        System.out.println("Detected non-terminal linear sugar moieties counter: " + tmpNonTerminalLinearSugarMoietiesCounter);
        System.out.println("Detected linear sugar moieties that are part of rings counter: " + tmpLinearSugarMoietiesInRingsCounter);
        System.out.println();
        System.out.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        System.out.println();
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        System.out.println("Basically a single sugar unit counter: " + tmpBasicallyASingleSugarUnitCounter);
        System.out.println("Basically a circular sugar counter: " + tmpBasicallyACircularSugarCounter);
        System.out.println("Basically a single circular sugar counter: " + tmpBasicallyASingleCircularSugarCounter);
        System.out.println("Basically a linear sugar counter: " + tmpBasicallyALinearSugarCounter);
        System.out.println("Basically a single linear sugar counter: " + tmpBasicallyASingleLinearSugarCounter);
        System.out.println();
        System.out.println("How many molecules have how many sugars: ");
        int tmpTotalOfSugarContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManySugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i));
            tmpTotalOfSugarContainingMolecules += tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("How many molecules have how many circular sugars: ");
        int tmpTotalOfCircularSugarContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i));
            tmpTotalOfCircularSugarContainingMolecules += tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("How many molecules have how many circular sugars attached via a glycosidic bond: ");
        int tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i));
            tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules += tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i);
        }
        System.out.println();
        System.out.println("How many molecules have how many linear sugars: ");
        int tmpTotalOfLinearSugarContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugarContainingMolecules += tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Size (= heavy atom count) frequency distribution of circular sugars: ");
        int tmpTotalOfCircularSugars = 0;
        for (int i : tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i));
            tmpTotalOfCircularSugars += tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Size (= heavy atom count) frequency distribution of linear sugars: ");
        int tmpTotalOfLinearSugars1 = 0;
        for (int i : tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugars1 += tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Size (= carbon atom count) frequency distribution of linear sugars: ");
        int tmpTotalOfLinearSugars2 = 0;
        for (int i : tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugars2 += tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Terminal Sugars counter: " + (tmpTerminalCircularSugarMoietiesCounter + tmpTerminalLinearSugarMoietiesCounter));
        System.out.println(tmpTerminalCircularSugarMoietiesCounter + " of these are circular");
        System.out.println(tmpTerminalLinearSugarMoietiesCounter + " of these are linear");
        System.out.println();
        System.out.println("Non-terminal Sugars counter: " + (tmpNonTerminalCircularSugarMoietiesCounter + tmpNonTerminalLinearSugarMoietiesCounter));
        System.out.println(tmpNonTerminalCircularSugarMoietiesCounter + " of these are circular");
        System.out.println(tmpNonTerminalLinearSugarMoietiesCounter + " of these are linear");
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atoms to atoms in ring ratios of circular sugars: ");
        List<String> tmpSortedKeySet = tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.keySet().stream().collect(Collectors.toList());
        Collections.sort(tmpSortedKeySet);
        for (String i : tmpSortedKeySet) {
            System.out.println(i + ": " + tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(i));
        }
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 5-membered circular sugars: ");
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(i));
        }
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 6-membered circular sugars: ");
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(i));
        }
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 7-membered circular sugars: ");
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(i));
        }
        System.out.println();
        System.out.println("Number of circular sugar moieties that had an unexpected ring size (should be zero!): " + tmpUnexpectedRingSizeCounter);
        System.out.println();
        System.out.println("Number of detected linear sugars in rings that got lost through the removal of circular " +
                "sugars counter: " + tmpLinearSugarsDetectedInCircularSugarsCounter);
        System.out.println();

        //Example print out of CNP list elements:
        //System.out.println("These molecules are basically sugars:");
        //System.out.println(Arrays.toString(tmpBasicallyASingleSugarUnitCNPs.toArray()));
        //System.out.println(tmpBasicallyASingleSugarUnitCNPs); //shorter alternative
        //</editor-fold>

        tmpCursor.close();

        //<editor-fold desc="Tests for consistency">
        Assert.assertEquals(tmpMoleculesCounter, tmpHasNoSugarsCounter + tmpHasAnyTypeOfSugarsCounter);
        Assert.assertEquals(tmpHasAnyTypeOfSugarsCounter, tmpHasOnlyCircularSugarsCounter
                + tmpHasOnlyLinearSugarsCounter
                + tmpHasCircularAndLinearSugarsCounter);
        Assert.assertEquals(tmpHasCircularSugarsCounter, tmpHasTerminalAndNonTerminalCircularSugarsCounter
                + tmpHasOnlyTerminalCircularSugarsCounter
                + tmpHasOnlyNonTerminalCircularSugarsCounter);
        Assert.assertTrue(tmpHasBoth >= 0);
        Assert.assertEquals(tmpHasLinearSugarsCounter, tmpHasTerminalAndNonTerminalLinearSugarsCounter
                + tmpHasOnlyTerminalLinearSugarsCounter
                + tmpHasOnlyNonTerminalLinearSugarsCounter);
        Assert.assertEquals(tmpHasAnyTypeOfSugarsCounter, tmpTotalOfSugarContainingMolecules);
        Assert.assertEquals(tmpHasCircularSugarsCounter, tmpTotalOfCircularSugarContainingMolecules);
        Assert.assertEquals(tmpHasGlycosidicBondCounter, tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules);
        Assert.assertEquals(tmpHasLinearSugarsCounter, tmpTotalOfLinearSugarContainingMolecules);
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTerminalCircularSugarMoietiesCounter
                + tmpNonTerminalCircularSugarMoietiesCounter);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTerminalLinearSugarMoietiesCounter
                + tmpNonTerminalLinearSugarMoietiesCounter);
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTotalOfCircularSugars);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTotalOfLinearSugars1);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTotalOfLinearSugars2);
        Assert.assertEquals(tmpBasicallyASingleSugarUnitCounter,
                tmpBasicallyASingleCircularSugarCounter + tmpBasicallyASingleLinearSugarCounter);
        //</editor-fold>
    }

    /**
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void linearSugarPatternsAppearanceTest() throws Exception {
        MongoClientSettings.Builder tmpBuilder = MongoClientSettings.builder();
        ServerAddress tmpAddress = new ServerAddress("localhost", 27017);
        tmpBuilder.applyToClusterSettings(builder -> builder.hosts(Collections.singletonList(tmpAddress)));
        MongoClientSettings tmpSettings = tmpBuilder.build();
        MongoClient tmpMongoClient = MongoClients.create(tmpSettings);
        String tmpCollectionName = "COCONUTmay";
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
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        // to also test their appearance
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        //Note: Here, additional molecules could be added to the list to also test them
        List<String> tmpLinearSugarsList = tmpSugarRemovalUtil.getLinearSugarStructuresList();
        List<List<Object>> tmpLinearSugarPatterns = new ArrayList<>(tmpLinearSugarsList.size());
        for (String tmpLinearSugarString : tmpLinearSugarsList) {
            List<Object> tmpList = new ArrayList<>(4);
            tmpList.add(0, tmpLinearSugarString);
            tmpList.add(1, DfPattern.findSubstructure(tmpSmiPar.parseSmiles(tmpLinearSugarString)));
            tmpList.add(2, 0);
            tmpLinearSugarPatterns.add(tmpList);
        }
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        int tmpMoleculeCounter = 0;
        int tmpContainsLinearSugarCounter = 0;
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
                tmpID = tmpCurrentDoc.getString("coconut_id");
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                boolean tmpMolHasAMatch = false;
                for (List<Object> tmpEntry : tmpLinearSugarPatterns) {
                    DfPattern tmpPattern = (DfPattern) tmpEntry.get(1);
                    if (tmpPattern.matches(tmpMolecule)) {
                        tmpEntry.set(2, (int)tmpEntry.get(2) + 1);
                        tmpMolHasAMatch = true;
                    }
                }
                if (tmpMolHasAMatch) {
                    tmpContainsLinearSugarCounter++;
                }
            } catch (Exception anException) {
                tmpLogger.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Molecules that have at least one match counter: " + tmpContainsLinearSugarCounter);
        for (List<Object> tmpEntry : tmpLinearSugarPatterns) {
           System.out.println((String)tmpEntry.get(0) + " " + (int)tmpEntry.get(2));
        }
        tmpCursor.close();
    }

    /**
     * TODO: Check whether the sugars or sugar-containing molecules in ZINC are in COCONUT
     */
    @Ignore
    @Test
    public void zincInspectionTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpZincSdfFile = new File(tmpClassLoader.getResource("zinc-all-for-sale.sdf").getFile());
        System.out.println(tmpZincSdfFile.getAbsolutePath());
        Logger tmpLogger = Logger.getLogger(SugarRemovalUtilityTest.class.getName());
        String tmpOutputFolderPath = (new File("SugarRemovalUtilityTest_Output")).getAbsolutePath() + File.separator
                + "zinc_inspection_test" + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
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
        Logger.getLogger("").setLevel(Level.WARNING);
        //Done for reproducibility
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        System.out.println(tmpSugarRemovalUtil.areOnlyCircularSugarsWithOGlycosidicBondDetected());
        System.out.println(tmpSugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeSetting());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeThresholdSetting());
        System.out.println(tmpSugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected());
        System.out.println(tmpSugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        System.out.println(tmpSugarRemovalUtil.areLinearSugarsInRingsDetected());
        System.out.println(tmpSugarRemovalUtil.arePropertiesAddedToSugarContainingMolecules());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMinSizeSetting());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting());
        System.out.println(tmpSugarRemovalUtil.areLinearAcidicSugarsDetected());
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpZincSdfFile), DefaultChemObjectBuilder.getInstance(), true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        int tmpContainsLinearSugarsCounter = 0;
        int tmpContainsCircularSugarsCounter = 0;
        int tmpBasicallyASugarCounter = 0;
        IAtomContainer tmpMolecule = null;
        String tmpID = null;
        while (tmpReader.hasNext()) {
            try {
                tmpMolecule = tmpReader.next();
                tmpID = tmpMolecule.getProperty("zinc_id");
                String tmpSmilesCode = tmpMolecule.getProperty("smiles");
                tmpMoleculeCounter++;
                tmpMolecule = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    /*tmpDepictionGenerator.depict(tmpSmiPar.parseSmiles(tmpSmilesCode)).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    /*tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                    tmpSugarContainingMoleculesCounter++;
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsCircularSugarsCounter++;
                    }
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        //tmpDepictionGenerator.depict(tmpSmiPar.parseSmiles(tmpSmilesCode)).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                        //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                        tmpContainsLinearSugarsCounter++;
                    }
                    /*if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true
                    && (boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        //tmpDepictionGenerator.depict(tmpSmiPar.parseSmiles(tmpSmilesCode)).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                        //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
                    }*/
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
        tmpReader.close();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Tests on specific molecules">
    /**
     * TODO
     * @throws Exception
     */
    @Test
    public void testCircularSugarRemoval() throws Exception {
        Map<String, String> tmpSmilesBeforeAndAfterDeglycosylationMap = new HashMap<>(3, 1.0f);
        tmpSmilesBeforeAndAfterDeglycosylationMap.put("CC(N)C(=O)NC(CCC(N)=O)C(=O)NOC1OC(O)C(O)C(O)C1O", //CHEMBL56258
                "O=C(N)CCC(NC(=O)C(N)C)C(=O)NO");
        tmpSmilesBeforeAndAfterDeglycosylationMap.put("CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O", //CHEMBL168422
                "O=C(OC1C(O)C(OC(CO)C1O)C=2C(O)=CC(O)=CC2CO)C=CC=CCC(O)C=CC=CCCCCC");
        tmpSmilesBeforeAndAfterDeglycosylationMap.put("OC1OC(O)C(O)C1OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1OC1C(O)C(O)C(O)OC(O)C1O", //own creation
                "OC1C(OCCCCCCCCCCC)OC(OCCCCCCCCCCCCCCCCC)C(O)C1O");
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
        IAtomContainer tmpOriginalMolecule;
        for (String tmpKey : tmpSmilesBeforeAndAfterDeglycosylationMap.keySet()) {
            tmpOriginalMolecule = tmpSmiPar.parseSmiles(tmpKey);
            Assert.assertTrue(tmpSugarRemovalUtil.hasCircularSugars(tmpOriginalMolecule));
            Assert.assertTrue(tmpSugarRemovalUtil.hasCircularAndOrLinearSugars(tmpOriginalMolecule));
            tmpOriginalMolecule = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, false);
            String tmpSmilesAfterDeglycosylation = tmpSmiGen.create(tmpOriginalMolecule);
            System.out.println(tmpSmilesAfterDeglycosylation);
            Assert.assertEquals(tmpSmilesBeforeAndAfterDeglycosylationMap.get(tmpKey), tmpSmilesAfterDeglycosylation);
        }
    }

    /**
     * Also, testing whether explicit hydrogens cause problems in the algorithm. They do not seem to do but the resulting
     * molecule will of course also have them
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("[H]OC1([H])C([H])(OC2=C3C(OC(=O)C4=C3C([H])([H])C([H])([H])C4([H])[H])=C([H])C(=C2[H])C([H])([H])[H])OC([H])(C([H])(O[H])C1([H])O[H])C([H])([H])O[H]"); //CNP0223287 in COCONUTmay
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //A simple example, the sugar has a glycosidic bond and is not terminal and therefore removed; The resulting
        // disconnected CH3OH is too small to keep and gets cleared away
        Assert.assertEquals("[H]C1=C(O)C2=C(OC(=O)C3=C2C([H])([H])C([H])([H])C3([H])[H])C([H])=C1C([H])([H])[H]", tmpSmilesCode);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4"); //CNP0000012 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The sugar ring is not terminal and should not be removed, so the molecule remains unchanged
        Assert.assertEquals("O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O"); //CNP0000006 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing is removed, the sugar is terminal because the two phosphate groups are big enough to keep
        Assert.assertEquals("O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O", tmpSmilesCode);
        tmpSugarRemovalUtil.setStructureToKeepModeThresholdSetting(6);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, one of the phosphate groups is removed because it has only 5 heavy atoms and therefore, the sugar is
        // no longer terminal and also removed
        Assert.assertEquals("O=P(O)(O)OC", tmpSmilesCode);
        tmpSugarRemovalUtil.setStructureToKeepModeThresholdSetting(7);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, both phosphate groups are removed because they are too small and nothing remains of the molecule
        Assert.assertEquals("", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        //back to default
        tmpSugarRemovalUtil.setStructureToKeepModeThresholdSetting(5);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, also non-terminal sugars are removed, which leaves two unconnected phosphate groups in this case
        Assert.assertEquals("O=P(O)(O)O.O=P(O)(O)OC", tmpSmilesCode);
        tmpSugarRemovalUtil.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(0.7);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, the sugar ring does not have enough oxygen atoms attached to be classified as a sugar and be removed
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC"); //CNP0000030 in COCONUTfebruary20
        Assert.assertFalse(tmpSugarRemovalUtil.hasLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //There is a structure at the center of the molecule that matches some of the linear sugar patterns (or did before)
        // but it should not be detected as a linear sugar and not be removed.
        Assert.assertEquals("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //now too, nothing should happen
        Assert.assertEquals("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // even with non-terminal sugars and linear sugars in rings removed, the structure should not change
        Assert.assertEquals("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC", tmpSmilesCode);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        // this molecule was chosen because before there was a bug detected here, a linear sugar candidate was detected
        // and a CH3 group removed
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=CCC12OC(OC1=O)C3(C)C(C)CCC3(C)C2"); //CNP0000023 in COCONUTfebruary20
        Assert.assertFalse(tmpSugarRemovalUtil.hasLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        // this molecule contains a macrocyle that is partly made up of sugars
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O"); //CNP0000082 in COCONUTfebruary20
        // since linear sugars in cycles should not be removed, this must be false
        Assert.assertFalse(tmpSugarRemovalUtil.hasLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here because the linear sugar in the macrocycle is not detected and it would not be terminal
        Assert.assertEquals("O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        Assert.assertTrue(tmpSugarRemovalUtil.hasLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here because the linear sugar in the macrocycle is detected but not terminal
        Assert.assertEquals("O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // now, the sugar in the macrocycle is removed and the molecule is disconnected
        Assert.assertEquals("OC=1C=CC=C(C1)C=2C=CC=C3OCC4C5=C(OC4C32)C(OC)=C(OC)C=C5CNC.NCCNC", tmpSmilesCode);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        // two structures connected by a circular sugar moiety
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=N[CH+]C5)C(C)C)C(O)C1O"); //CNP0000233 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here because the circular sugar is not terminal
        Assert.assertEquals("O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=N[CH+]C5)C(C)C)C(O)C1O", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed because only linear sugars are removed and there are none!
        Assert.assertEquals("O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=N[CH+]C5)C(C)C)C(O)C1O", tmpSmilesCode);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, the circular sugar is removed
        Assert.assertEquals("O=C1C=2C=CC(O)=CC2[O+]=C[C-]1C3=CC=C(O)C=C3.N1=CC(=C[CH+]1)CC(CNC)C(C)C", tmpSmilesCode);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC(C2=COC=C2)CC3(C)C1CCC4(C)C3C5OC(=O)C4(O)C=C5"); //CNP0000304 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there might be a match for the linear sugar patterns
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC2CC3(OC4(O)C(CC5(OC45C(=O)OC)CCCCCCCCCCCCCCCC)C2(O3)C1)CCCCCCCCCCCCCCCC"); //CNP0000445 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is might be match for the linear sugar patterns
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        //Interesting case because there is a macrocycle containing a sugar cycle that is not isolated
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO"); //CNP0000509 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here although there is a match for the linear sugar patterns in the said non-isolated sugar cycle in the macrocycle
        Assert.assertEquals("O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, the sugar cycle in the macrocyle is removed because it matches the linear sugar patterns
        Assert.assertEquals("O=C1C2=CC=CC(=C2CN1CC(=O)C3=C(O)C=4C=C(C=CCC5=CNC=6C=CC(=CC65)CNC)C(OC)=C7C4C=8C(=CC9CCCC9C38)CC7)CC#C", tmpSmilesCode);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(OC1OC(CO)C(O)C(O)C1O)(C)CC(=O)OCC=CC2=CC(OC)=C(OC3OC(CO)C(O)C(O)C3O)C(OC)=C2"); //CNP0000920 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The two circular sugar moieties and one connected linear sugar (a sugar acid) are removed (all terminal)
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1=CC=C(OC2OC(CO)C(O)C(O)C2O)C=C1"); //CNP0001189 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //One circular and one linear sugar (a sugar acid) moiety are removed
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(COC(=O)C(CC(=O)O)C(OCC3C(=C)CCC4C(C)(C)CCCC34C)C(=O)OC)CCCC12C)C(=O)OC"); //CNP0001863 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //One linear sugar moiety is removed (a sugar acid), the other one is not because it is terminal; no other changes are done to the molecule
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2OC%10OC(CO)C(O)C(O)C%10O)CO)C(O)C1O"); //CNP0002871 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println();
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("OC1CC2CCC3C(CCC4(C)C3CC5OC6(OCC(C)CC6)C(C)C54)C2(C)CC1O", tmpSmilesCode);
        //for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMoleculesAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC4C5=CCC6C7(C)CCC(O)C(C(=O)OC8OC(CO)C(O)C(O)C8O)(C)C7CCC6(C)C5(C)CC3)C(O)C(OC9OC(CO)C(O)C(O)C9O)C2O)C(OC%10OC(CO)C(O)C(O)C%10O)C(O)C1O"); //CNP0005247 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println();
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(O)C1(C)C(O)CCC2(C)C1CCC3(C)C2CC=C4C5CC(C)(C)CCC5(C(=O)O)CCC43C", tmpSmilesCode);
        //for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMoleculesAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OC1C(O)C(OC2C3=C(O)C(=CC=C3OC2C(=C)CO)C(=O)C)OC(CO)C1O"); //CNP0032326 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println();
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(C1=CC=C2OC(C(=C)CO)C(O)C2=C1O)C", tmpSmilesCode);
        //for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMoleculesAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"); //CNP0031401 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println();
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C1C=C(OC=2C=C(O)C=C(O)C12)C=3C=CC(O)=C(O)C3", tmpSmilesCode);
        //for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMoleculesAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2OC3OCC(O)C(OC4OCC(O)C(O)C4O)C3O)OC5C(OC(=O)C67CCC(C)(C)CC7C8=CCC9C%10(C)CC(O)C(OC%11OC(CO)C(O)C(O)C%11O)C(C(=O)O)(C)C%10CCC9(C)C8(CO)CC6)OC(C)C(OC(=O)C=CC%12=CC(OC)=C(OC)C(OC)=C%12)C5OC%13OC(C)C(O)C(O)C%13O)C1O"); //CNP0028122 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println();
        System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before), only one non-terminal sugar remains
        Assert.assertEquals("O=C(OC1C(O)C(O)C(OC(=O)C23CCC(C)(C)CC3C4=CCC5C6(C)CC(O)C(O)C(C(=O)O)(C)C6CCC5(C)C4(CO)CC2)OC1C)C=CC7=CC(OC)=C(OC)C(OC)=C7", tmpSmilesCode);
        //for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMoleculesAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C1OC(O)C(O)C(O)C1O"); //CNP0002919 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //This molecule is a sugar and should be completely removed
        Assert.assertEquals("", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Even with detection of glycosidic bonds, this sugar ring should be removed because there is no other structure
        // it can bind to via a glycosidic bond
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(C)CCCC12C)C(=O)O"); //CNP0003329 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The linear sugar moiety (a sugar di-acid) is removed
        Assert.assertEquals("C=C1CCC2C(C)(C)CCCC2(C)C1C", tmpSmilesCode);
        //another example CNP0031156 in COCONUTfebruary20
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        // interesting is the peroxide bond here
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C=C(OC=2C1=CC3=C(OC(C)(C)C(OOCC(O)C(O)C(O)C(O)CO)C3)C2[N+]=4C=C5N=CC=C5C4CC)C"); //CNP0007654 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4"); //CNP0032817 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The sugar moiety is NOT connected to the core structure via a glycosidic bond, so it is not removed
        Assert.assertEquals("O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that this setting is changed, the sugar moiety is removed
        //note: chemically, the carboxy group should be part of the sugar, not of the core
        Assert.assertEquals("O=COC=1C=C2OC(=CC(=O)C2=C(O)C1O)C=3C=CC(O)=CC3", tmpSmilesCode);

        //another examples for the same thing:
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"); //CNP0031401 in COCONUTfebruary20
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The circular sugar moiety is NOT connected to the core structure via a glycosidic bond, so it is not removed
        //also, the removal of linear sugars leaves the circular sugar untouched
        Assert.assertEquals("O=C1C=C(OC=2C1=C(O)C=C(O)C2C3OC(CO)C(O)C(O)C3O)C=4C=CC(O)=C(O)C4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC"); //CNP0000913 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Only the aliphatic chain remains (in its full length)
        Assert.assertEquals("CCCCCCCCCCCCCC", tmpSmilesCode);
    }


    @Test
    public void specificTest20() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=CC1(C)C(OC2OC(C(=O)O)C(O)C(OC3OCC(O)C(O)C3O)C2OC4OC(CO)C(O)C(O)C4O)CCC5(C)C6CC=C7C8CC(C)(C)CCC8(C(=O)OC9OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(OC%10OC(CO)C(O)C%10O)C(C)CC)C(C)CC)C(O)C9OC%11OC(C)C(OC%12OCC(O)C(O)C%12O)C(O)C%11O)C(O)CC7(C)C6(C)CCC15"); //CNP0000306 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //Only the core structure remains
        Assert.assertEquals("O=CC1(C)C(O)CCC2(C)C1CCC3(C)C2CC=C4C5CC(C)(C)CCC5(C(=O)O)C(O)CC43C", tmpSmilesCode);
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMoleculesAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));
    }

    @Test
    public void specificTest21() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(OCC)CCC1=CC=2C=COC2C=3OCCNCC4(SSCC5(OC(OC31)C(O)C(O)C5O)CO)CCCC4"); //CNP0001379 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //the sugar ring within the macrocycle can be matched by the linear sugar patterns but should not be removed
        Assert.assertEquals("O=C(OCC)CCC1=CC=2C=COC2C=3OCCNCC4(SSCC5(OC(OC31)C(O)C(O)C5O)CO)CCCC4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, the sugar (circular) sugar moiety within the macrocyle is removed
        Assert.assertEquals("O=C(OCC)CCC=1C=C(OCCNCC2(SSC)CCCC2)C=3OC=CC3C1", tmpSmilesCode);
    }

    @Test
    public void specificTest22() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("CCC(CC1CCCC2(O1)CC3C(C(O2)CC4(C(CC(O4)(C)C)C=CCCCCCC(C(C(C(C(C(C(C=CC(=O)O3)(C)O)O)C)O)OC5CCC(C(O5)C)N(C)C)O)(C)O)O)C)O"); //Ossamycin
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing gets removed
        Assert.assertEquals("O=C1OC2CC3(OC(CCC3)CC(O)CC)OC(CC4(O)OC(C)(C)CC4C=CCCCCCC(O)(C)C(O)C(OC5OC(C)C(N(C)C)CC5)C(O)C(C)C(O)C(O)(C=C1)C)C2C", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Again, nothing is removed because the linear sugar is non-terminal AND contained in a macrocycle
        Assert.assertEquals("O=C1OC2CC3(OC(CCC3)CC(O)CC)OC(CC4(O)OC(C)(C)CC4C=CCCCCCC(O)(C)C(O)C(OC5OC(C)C(N(C)C)CC5)C(O)C(C)C(O)C(O)(C=C1)C)C2C", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, the linear sugar gets removed
        Assert.assertEquals("O=C(OC1CC2(OC(CCC2)CC(O)CC)OC(CC3(O)OC(C)(C)CC3C=CCCCCC)C1C)C=CC(O)(C)C(O)CC.O1CCCC(N(C)C)C1C", tmpSmilesCode);
    }

    /**
     * Interesting, because the linear sugar patterns do not match the entire sugar ring because of the amine. Therefore,
     * a part of the ring ends up in the linear sugar candidates (false-positive).
     */
    @Test
    public void specificTest23() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C"); //CNP0000225 in COCOCNUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed because this molecule has only circular sugars (was a problem before)
        Assert.assertEquals("O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C", tmpSmilesCode);
    }

    /**
     *
     */
    @Test
    public void specificTest24() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        // a molecule containing 2 circular sugars and 2 linear sugars (sugar acids), in both cases 1 terminal and 1 non-terminal
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(C)(O)CC(=O)OCc1ccc(cc1)OC1C(CO)OC(OC(C(=O)OCc2ccc(OC3OC(CO)CC(O)C3O)cc2)C(O)(CC(C)C)C(=O)OC2CCc3cc4cc(O)c(C)c(O)c4c(O)c3C2=O)C(O)C1O"); //imaginary molecule created for demonstration purposes
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // the terminal linear and circular sugar moieties are removed
        Assert.assertEquals("O=C(OCC1=CC=C(O)C=C1)C(OC2OC(CO)C(OC3=CC=C(C=C3)C)C(O)C2O)C(O)(C(=O)OC4C(=O)C5=C(O)C6=C(O)C(=C(O)C=C6C=C5CC4)C)CC(C)C", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // all 4 sugar moieties are removed, the two non-terminal sugar moieties included
        Assert.assertEquals("O=C1C2=C(O)C3=C(O)C(=C(O)C=C3C=C2CCC1)C.OC1=CC=C(C=C1)C.OC1=CC=C(C=C1)C", tmpSmilesCode);
        // back to default setting
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // only the terminal linear sugar is removed
        Assert.assertEquals("O=C(OCC1=CC=C(OC2OC(CO)CC(O)C2O)C=C1)C(OC3OC(CO)C(OC4=CC=C(C=C4)C)C(O)C3O)C(O)(C(=O)OC5C(=O)C6=C(O)C7=C(O)C(=C(O)C=C7C=C6CC5)C)CC(C)C", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // the two linear sugar moieties are removed, disconnecting the molecule
        Assert.assertEquals("O=C1C2=C(O)C3=C(O)C(=C(O)C=C3C=C2CCC1)C.OCC1OC(OC2=CC=C(C=C2)C)C(O)C(O)C1.OCC1OC(O)C(O)C(O)C1OC2=CC=C(C=C2)C", tmpSmilesCode);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // the two circular sugar moieties are removed, disconnecting the molecule
        Assert.assertEquals("O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(O)C=C1.O=C(OCC1=CC=C(O)C=C1)C(O)C(O)(C(=O)OC2C(=O)C3=C(O)C4=C(O)C(=C(O)C=C4C=C3CC2)C)CC(C)C", tmpSmilesCode);
        // back to default setting
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // only the terminal circular sugar moiety is removed
        Assert.assertEquals("O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(OC2C(O)C(O)C(OC(C(=O)OCC3=CC=C(O)C=C3)C(O)(C(=O)OC4C(=O)C5=C(O)C6=C(O)C(=C(O)C=C6C=C5CC4)C)CC(C)C)OC2CO)C=C1", tmpSmilesCode);
    }

    /**
     * TODO
     */
    @Test
    public void specificTest25() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("OC1(OCCC21OC3(O)CCOC3(O2)C)C"); //CNP0017323 in COCONUTmay
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed because the sugar is a spiro ring (before, these were not filtered)
        Assert.assertEquals("OC1(OCCC21OC3(O)CCOC3(O2)C)C", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that the spiro rings are detected, the sugar ring is removed; BUT the adjacent ring is left intact!
        Assert.assertEquals("OC12OCOC2(OCC1)C", tmpSmilesCode);
        //back to default for the other tests
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(false);

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("OCC1OC2(OC3C(O)C(OC3(OC2)CO)CO)C(O)C1O"); //CNP0309652 in COCONUTmay
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed because the sugar is a spiro ring (before, these were not filtered)
        Assert.assertEquals("OCC1OC2(OC3C(O)C(OC3(OC2)CO)CO)C(O)C1O", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that the spiro rings are detected, the sugar ring is removed; BUT the adjacent ring is left intact!
        Assert.assertEquals("OCC1OC2(OCCOC2C1O)CO", tmpSmilesCode);
        //back to default for the other tests
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(false);

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("OCC1OC2(OCC3OC(OC2)(CO)C(O)C3O)C(O)C1O"); //CNP0344838 in COCONUTmay
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed because the sugar is a spiro ring (before, these were not filtered)
        Assert.assertEquals("OCC1OC2(OCC3OC(OC2)(CO)C(O)C3O)C(O)C1O", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that the spiro rings are detected, the sugar ring is removed; BUT the adjacent ring is left intact!
        Assert.assertEquals("OCC12OCCOCC(O1)C(O)C2O", tmpSmilesCode);
        //back to default for the other tests (if some are added later...)
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(false);
    }

    /**
     * TODO
     */
    @Test
    public void specificTest26() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();

        //A molecule consisting entirely of one circular and one linear sugar
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O"); //CNP0121254 in COCONUTmay
        Assert.assertEquals(2, tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpOriginalMolecule));
        Assert.assertEquals(1, tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpOriginalMolecule));
        Assert.assertEquals(1, tmpSugarRemovalUtil.getNumberOfLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Completely removed
        Assert.assertEquals("", tmpSmilesCode);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Only the linear sugar remains
        Assert.assertEquals("OCC(O)C(O)C(O)CO", tmpSmilesCode);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Only the circular sugar remains
        Assert.assertEquals("OCC1OCC(O)C(O)C1O", tmpSmilesCode);

        //A molecule consisting entirely of one circular and one linear sugar
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("OCC(O)C(O)C(O)C(O)C(O)C1OC(O)C(O)C(O)C1N"); //CNP0185588 in COCONUTmay
        Assert.assertEquals(2, tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpOriginalMolecule));
        Assert.assertEquals(1, tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpOriginalMolecule));
        Assert.assertEquals(1, tmpSugarRemovalUtil.getNumberOfLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Completely removed
        Assert.assertEquals("", tmpSmilesCode);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Only the linear sugar remains
        Assert.assertEquals("OCC(O)C(O)C(O)C(O)CO", tmpSmilesCode);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Only the circular sugar remains
        Assert.assertEquals("OC1OCC(N)C(O)C1O", tmpSmilesCode);
    }

    /**
     * TODO tidy up!
     */
    @Test
    public void specificTest27() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpSugarRemovalUtil.setLinearSugarCandidateMaxSizeSetting(36);

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("C(C1C2C(C(C(O1)OC3C(OC(C(C3O)O)OC4C(OC(C(C4O)O)OC5C(OC(C(C5O)O)OC6C(OC(C(C6O)O)OC7C(OC(O2)C(C7O)O)CO)CO)CO)CO)CO)O)O)O"); //alpha-cyclodextrin
        RingSearch tmpRingSearch = new RingSearch(tmpOriginalMolecule);
        System.out.println(tmpRingSearch.isolatedRingFragments().size());
        List<IAtomContainer> tmpCandidateList = tmpSugarRemovalUtil.splitEtherEsterAndPeroxideBonds(tmpSugarRemovalUtil.combineOverlappingCandidates(tmpSugarRemovalUtil.detectLinearSugarCandidatesByPatternMatching(tmpOriginalMolecule)));
        tmpSugarRemovalUtil.removeAtomsOfCircularSugarsFromCandidates(tmpCandidateList, tmpOriginalMolecule);
        for (IAtomContainer tmpCandidate : tmpCandidateList) {
            System.out.println(tmpSmiGen.create(tmpCandidate));
        }
        System.out.println(tmpSugarRemovalUtil.hasCircularAndOrLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //
        Assert.assertEquals("", tmpSmilesCode);
    }
    //</editor-fold>

    /**
     *
     */
    @Test
    public void classPropertiesTest() throws Exception {
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true); //default would be false
        tmpSugarRemovalUtil.addLinearSugarToStructureList("CO");
        List<String> tmpLinearSugarsList = tmpSugarRemovalUtil.getLinearSugarStructuresList();
        for (String tmpSmiles : tmpLinearSugarsList) {
            System.out.println(tmpSmiles);
        }
        System.out.println("-------------------");
        tmpSugarRemovalUtil.addCircularSugarToStructureList("O1CCCCCCC1");
        List<String> tmpCircularSugarsList = tmpSugarRemovalUtil.getCircularSugarStructuresList();
        for (String tmpSmiles : tmpCircularSugarsList) {
            System.out.println(tmpSmiles);
        }
        System.out.println("-------------------");
        System.out.println(tmpSugarRemovalUtil.areOnlyCircularSugarsWithOGlycosidicBondDetected());
        System.out.println(tmpSugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeSetting());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeThresholdSetting());
        System.out.println(tmpSugarRemovalUtil.areOnlyCircularSugarsWithEnoughExocyclicOxygenAtomsDetected());
        System.out.println(tmpSugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        System.out.println(tmpSugarRemovalUtil.areLinearSugarsInRingsDetected());
        System.out.println(tmpSugarRemovalUtil.arePropertiesAddedToSugarContainingMolecules());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMinSizeSetting());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMaxSizeSetting());
        System.out.println(tmpSugarRemovalUtil.areLinearAcidicSugarsDetected());
        System.out.println(tmpSugarRemovalUtil.areSpiroRingsDetectedAsCircularSugars());
    }

    //<editor-fold desc="Tests for protected routines">
    /**
     * TODO: make this a real test with example structures to split
     */
    @Test
    public void testEtherEsterPeroxideSplitting() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer target1 = tmpSmiPar.parseSmiles("CCOCC");
        IAtomContainer target2 = tmpSmiPar.parseSmiles("CC(=O)OCC");
        IAtomContainer target3 = tmpSmiPar.parseSmiles("CCOOCC");
        List<IAtomContainer> tmpCandidates = new ArrayList<>(3);
        tmpCandidates.add(target1);
        tmpCandidates.add(target2);
        tmpCandidates.add(target3);
        for (IAtomContainer tmpCandidate : this.splitEtherEsterAndPeroxideBonds(tmpCandidates)) {
            System.out.println(tmpSmiGen.create(tmpCandidate));
        }
    }
    //</editor-fold>

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

    @Test
    public void experiments() throws Exception {
        IAtomContainer tmpMol = new AtomContainer();
        IAtom tmpAtom = new Atom();
        tmpMol.addAtom(tmpAtom);
        tmpMol.addAtom(tmpAtom);
        System.out.println(tmpMol.getAtomCount());
        NumberFormat df = NumberFormat.getInstance(Locale.US);
        df.setMaximumFractionDigits(1);
        df.setRoundingMode(RoundingMode.DOWN);
        System.out.println(df.format(4.99));
    }

    @Test
    public void smartsTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer target = tmpSmiPar.parseSmiles("CCOCC");
        SmartsPattern tmpEtherPattern = SmartsPattern.create("[CD2]-[OX2]-[CD2]");
        Mappings tmpMappings = tmpEtherPattern.matchAll(target).uniqueAtoms();
        for (IAtomContainer tmpMap : tmpMappings.toSubstructures()) {
            System.out.println(tmpSmiGen.create(tmpMap));
            IAtom tmpCarbon1 = null;
            IAtom tmpCarbon2 = null;
            IAtom tmpOxygen = null;
            for (IAtom tmpAtom : tmpMap.atoms()) {
                if (tmpAtom.getSymbol().equals("O")) {
                    tmpOxygen = tmpAtom;
                } else if (tmpAtom.getSymbol().equals("C") && Objects.isNull(tmpCarbon1)) {
                    tmpCarbon1 = tmpAtom;
                } else {
                    tmpCarbon2 = tmpAtom;
                }
            }
            target.removeBond(tmpOxygen, tmpCarbon2);
            tmpOxygen.setImplicitHydrogenCount(1);
            IAtom tmpNewOxygen = new Atom("O");
            tmpNewOxygen.setImplicitHydrogenCount(1);
            target.addAtom(tmpNewOxygen);
            IBond tmpNewBond = new Bond(tmpNewOxygen, tmpCarbon2, IBond.Order.SINGLE);
            target.addBond(tmpNewBond);
        }
        IAtomContainerSet tmpComponents = ConnectivityChecker.partitionIntoMolecules(target);
        for (IAtomContainer tmpComponent : tmpComponents.atomContainers()) {
            System.out.println(tmpSmiGen.create(tmpComponent));
        }
        target = tmpSmiPar.parseSmiles("CC(=O)OCC");
        SmartsPattern tmpEsterPattern = SmartsPattern.create("[CD3](=[OX1])-[OX2]-[CD2]");
        tmpMappings = tmpEsterPattern.matchAll(target).uniqueAtoms();
        for (IAtomContainer tmpMap : tmpMappings.toSubstructures()) {
            System.out.println(tmpSmiGen.create(tmpMap));
        }

    }

    @Test
    public void testIsomorphismDependencyOnValences() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMissingValenceRing = tmpSmiPar.parseSmiles("O=C(O)C1OC([O])[C](O)C(O)C1O");
        RingSearch tmpRingSearch = new RingSearch(tmpMissingValenceRing);
        IAtomContainer tmpReferenceRing = tmpSmiPar.parseSmiles("C1CCOCC1");
        UniversalIsomorphismTester tmpUnivIsoTester = new UniversalIsomorphismTester();
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        System.out.println(tmpSmiGen.create(tmpRingSearch.isolatedRingFragments().get(0)));
        boolean tmpAreIsomorph = tmpUnivIsoTester.isIsomorph(tmpRingSearch.isolatedRingFragments().get(0), tmpReferenceRing);
        System.out.println(tmpAreIsomorph);
    }

    @Test
    public void testEtherEsterPeroxideSplitting2() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpCandidate = tmpSmiPar.parseSmiles("O=C(O)C[C](O)CC(=O)OC1COC(OC2C(O)C(OC(CO)C2O)OC3C(O)C(O)C([O])OC3CO)C(O)C1O");
        SmartsPattern tmpEtherPattern = SmartsPattern.create("[C]-[OX2;!R]-[C]");
        Mappings tmpMapping = tmpEtherPattern.matchAll(tmpCandidate).uniqueAtoms();
        System.out.println(tmpMapping.count());
        List<IAtomContainer> tmpList = new ArrayList<>(1);
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        tmpList.add(tmpCandidate);
        for (IAtomContainer tmpNewCandidate : this.splitEtherEsterAndPeroxideBonds(tmpList)) {
            System.out.println(tmpSmiGen.create(tmpNewCandidate));
        }
    }

    /**
     * Molecule from specificTest10
     * @throws Exception
     */
    @Test
    public void depictionTest() throws Exception {
        String tmpOutputFolderPath = (new File("SugarRemovalUtilityTest_Output")).getAbsolutePath() + File.separator
                + "easy_depiction_test" + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
        //Interesting case because there is a macrocycle containing a sugar cycle that is not isolated
        //tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO"); //CNP0000509 in COCONUTfebruary20
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1"); //CNP0256712 in COCONUTmay
        //tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=CC(O)C(O)C(O)C(O)COC(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O"); //CNP0140380 in COCONUTmay
        tmpDepictionGenerator.withSize(2000, 2000)
                .withFillToFit()
                .depict(tmpOriginalMolecule)
                .writeTo(tmpOutputFolderPath + File.separator + "Test_original_molecule.png");
        List<IAtomContainer> tmpCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpOriginalMolecule);
        List<IAtomContainer> tmpToHighlight = new ArrayList<>(tmpCandidates.size());
        for (int i = 0; i < tmpCandidates.size(); i++) {
            IAtomContainer tmpCandidate = tmpCandidates.get(i);
            tmpToHighlight.add(tmpCandidate);
        }
        tmpDepictionGenerator.withHighlight(tmpToHighlight, Color.BLUE)
                .withSize(2000, 2000)
                .withFillToFit()
                .depict(tmpOriginalMolecule)
                .writeTo(tmpOutputFolderPath + File.separator + "Test" + ".png");
    }

    //<editor-fold desc="Protected methods">
    /**
     * TODO
     * @return
     */
    protected SugarRemovalUtility getSugarRemovalUtilityV0100DefaultSettings() {
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
        tmpSugarRemovalUtil.setStructureToKeepModeSetting(StructureToKeepModeOption.HEAVY_ATOM_COUNT);
        tmpSugarRemovalUtil.setStructureToKeepModeThresholdSetting(5);
        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(true);
        tmpSugarRemovalUtil.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(0.5);
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
        tmpSugarRemovalUtil.setLinearSugarCandidateMinSizeSetting(4);
        tmpSugarRemovalUtil.setLinearSugarCandidateMaxSizeSetting(7);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(false);
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(false);
        return tmpSugarRemovalUtil;
    }
    //</editor-fold>
}
