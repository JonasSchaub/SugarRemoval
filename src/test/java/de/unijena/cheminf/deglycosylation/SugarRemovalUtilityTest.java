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

import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

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
        System.out.println(tmpSugarRemovalUtil.isGlycosidicBondDetected());
        System.out.println(tmpSugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println(tmpSugarRemovalUtil.getStructuresToKeepMode());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeThreshold());
        System.out.println(tmpSugarRemovalUtil.isNrOfAttachedOxygensIncluded());
        System.out.println(tmpSugarRemovalUtil.getAttachedOxygensToAtomsInRingRatioThreshold());
        System.out.println(tmpSugarRemovalUtil.areLinearSugarsInRingsRemoved());
        System.out.println(tmpSugarRemovalUtil.arePropertiesOfSugarContainingMoleculesSet());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMinSize());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMaxSize());
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
                tmpMolecule = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    //tmpDepictionGenerator.depict(tmpSmiPar.parseSmiles(tmpSmilesCode)).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    //tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_1.png");
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
        tmpCursor.close();
    }

    /**
     *
     */
    @Ignore
    @Test
    public void CoconutStatsTest() throws Exception {
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
        //Done for reproducibility
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        List<String> tmpSugarContainingMoleculesCNPs = new ArrayList<>(50000);
        int tmpContainsCircularSugarsCounter = 0;
        List<String> tmpContainsCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasTerminalRingsCounter = 0;
        List<String> tmpHasTerminalRingsCNPS = new ArrayList(50000);
        int tmpHasNonTerminalRingsCounter = 0;
        List<String> tmpHasNonTerminalRingsCNPS = new ArrayList(50000);
        int tmpHasOnlyTerminalRingsCounter = 0;
        List<String> tmpHasOnlyTerminalRingsCNPS = new ArrayList(50000);
        int tmpHasOnlyNonTerminalRingsCounter = 0;
        List<String> tmpHasOnlyNonTerminalRingsCNPS = new ArrayList(50000);
        int tmpHasOnlyCircularSugarsCounter = 0;
        List<String> tmpHasOnlyCircularSugarsCNPS = new ArrayList(50000);
        int tmpHasGlycosidicBondCounter = 0;
        List<String> tmpHasGlycosidicBondCNPS = new ArrayList(50000);
        int tmpHasGlycosidicBondOnTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnTerminalSugarCNPS = new ArrayList(50000);
        int tmpHasGlycosidicBondOnNonTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnNonTerminalSugarCNPS = new ArrayList(50000);
        int tmpContainsLinearSugarsCounter = 0; //total number of molecules with linear sugars
        List<String> tmpContainsLinearSugarsCNPS = new ArrayList(3000);
        int tmpHasTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalLinearSugarsCNPS = new ArrayList(3000);
        int tmpHasNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasNonTerminalLinearSugarsCNPS = new ArrayList(3000);
        int tmpHasOnlyTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalLinearSugarsCNPS = new ArrayList(3000);
        int tmpHasOnlyNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalLinearSugarsCNPS = new ArrayList(3000);
        int tmpHasOnlyLinearSugarsCounter = 0;
        List<String> tmpHasOnlyLinearSugarsCNPS = new ArrayList(3000);
        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList(2000);
        int tmpBasicallyASugarCounter = 0; //number of molecules that are basically sugars (circular and linear combined)
        List<String> tmpBasicallyASugarCNPS = new ArrayList(2000);
        int tmpBasicallyACircularSugarCounter = 0;
        List<String> tmpBasicallyACircularSugarCNPS = new ArrayList(2000);
        int tmpBasicallyALinearSugarCounter = 0;
        List<String> tmpBasicallyALinearSugarCNPS = new ArrayList(2000);
        int tmpBasicallyASingleSugarUnitCounter = 0;
        List<String> tmpBasicallyASingleSugarUnitCNPs = new ArrayList(1000);
        int tmpBasicallyASingleCircularSugarCounter = 0;
        List<String> tmpBasicallyASingleCircularSugarCNPs = new ArrayList(1000);
        int tmpBasicallyASingleLinearSugarCounter = 0;
        List<String> tmpBasicallyASingleLinearSugarCNPs = new ArrayList(1000);
        int tmpCircularSugarsCounter = 0;
        int tmpCircularSugarsWithGlycosidicBondCounter = 0; //the rest have no glycosidic bond
        int tmpTerminalCircularSugarsWithGlycosidicBondCounter = 0; //the rest are non-terminal
        //TODO!!!
        HashMap<Integer, Integer> tmpSizesOfCircularSugarsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpSizesOfLinearSugarsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManySugarsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyLinearSugarsMap = new HashMap(10, 0.9f);
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
                tmpID = tmpCurrentDoc.getString("coconut_id");
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //using default settings (of version 0.1.0.0) where nothing else is specified
                boolean tmpHasAnyTypeOfSugar = tmpSugarRemovalUtil.hasCircularAndOrLinearSugars(tmpMolecule);
                // note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                // note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    tmpSugarContainingMoleculesCounter++;
                    tmpSugarContainingMoleculesCNPs.add(tmpID);
                    int tmpNumberOfCircularAndLinearSugarMoieties = tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpMolecule);
                    Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManySugarsMap.get(tmpNumberOfCircularAndLinearSugarMoieties);
                    if (Objects.isNull(tmpCurrentListValue)) {
                        tmpHowManyMoleculesHaveHowManySugarsMap.put(tmpNumberOfCircularAndLinearSugarMoieties, 1);
                    } else {
                        tmpHowManyMoleculesHaveHowManySugarsMap.put(tmpNumberOfCircularAndLinearSugarMoieties, tmpCurrentListValue + 1);
                    }
                    if (tmpHasAnyCircularSugar) {
                        tmpContainsCircularSugarsCounter++;
                        tmpContainsCircularSugarsCNPs.add(tmpID);
                        if (!tmpHasAnyLinearSugar) {
                            tmpHasOnlyCircularSugarsCounter++;
                            tmpHasOnlyCircularSugarsCNPS.add(tmpID);
                        }
                        //terminal and non-terminal, having a glycosidic bond or not (see default settings)
                        List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule);
                        int tmpAllCircularMoietiesCounter;
                        int tmpTerminalCounter;
                        int tmpNonTerminalCounter;
                        int tmpGlycosidicBondCounter;
                        int tmpGlycosidicBondAndTerminalCounter;
                        int tmpGlycosidicBondAndNonTerminalCounter;
                        tmpAllCircularMoietiesCounter = tmpCircularSugarCandidatesList.size();
                        Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyCircularSugarsMap.get(tmpAllCircularMoietiesCounter);
                        if (Objects.isNull(tmpCurrentValue)) {
                            tmpHowManyMoleculesHaveHowManyCircularSugarsMap.put(tmpAllCircularMoietiesCounter, 1);
                        } else {
                            tmpHowManyMoleculesHaveHowManyCircularSugarsMap.put(tmpAllCircularMoietiesCounter, tmpCurrentValue + 1);
                        }
                        //the first is the overall counter, the second one is specific for this molecule
                        tmpCircularSugarsCounter += tmpAllCircularMoietiesCounter;
                        //note: circular moieties that become terminal after removal of a linear moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalCircularMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        // -1 for the deglycosylated core at the beginning of the list
                        tmpTerminalCounter = tmpRemovedTerminalCircularMoieties.size() - 1 ;
                        tmpNonTerminalCounter = tmpAllCircularMoietiesCounter - tmpTerminalCounter;
                        Assert.assertTrue(tmpNonTerminalCounter >= 0);
                        // leaving default! Now, only circular sugars having glycosidic bonds are in the candidates and removed moieties
                        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
                        List<IAtomContainer> tmpCircularSugarCandidatesWithGlycosidicBondsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule);
                        tmpGlycosidicBondCounter = tmpCircularSugarCandidatesWithGlycosidicBondsList.size();
                        List<IAtomContainer> tmpRemovedTerminalCircularMoietiesWithGlycosidicBond = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        tmpGlycosidicBondAndTerminalCounter = tmpRemovedTerminalCircularMoietiesWithGlycosidicBond.size() - 1;
                        tmpGlycosidicBondAndNonTerminalCounter = tmpGlycosidicBondCounter - tmpGlycosidicBondAndTerminalCounter;
                        Assert.assertTrue(tmpGlycosidicBondAndNonTerminalCounter >= 0);
                        // back to default!
                        tmpSugarRemovalUtil.setDetectGlycosidicBond(false);
                        if (tmpTerminalCounter > 0) {
                            tmpHasTerminalRingsCounter++;
                            tmpHasTerminalRingsCNPS.add(tmpID);
                            if (tmpNonTerminalCounter == 0) {
                                tmpHasOnlyTerminalRingsCounter++;
                                tmpHasOnlyTerminalRingsCNPS.add(tmpID);
                            }
                        }
                        if (tmpNonTerminalCounter > 0) {
                            tmpHasNonTerminalRingsCounter++;
                            tmpHasNonTerminalRingsCNPS.add(tmpID);
                            if (tmpTerminalCounter == 0) {
                                tmpHasOnlyNonTerminalRingsCounter++;
                                tmpHasOnlyNonTerminalRingsCNPS.add(tmpID);
                            }
                        }
                        if (tmpGlycosidicBondCounter > 0) {
                            tmpHasGlycosidicBondCounter++;
                            tmpHasGlycosidicBondCNPS.add(tmpID);
                            //the first is the overall counter, the second one is specific for this molecule
                            tmpCircularSugarsWithGlycosidicBondCounter += tmpGlycosidicBondCounter;
                            if (tmpGlycosidicBondAndTerminalCounter > 0) {
                                tmpHasGlycosidicBondOnTerminalSugarCounter++;
                                tmpHasGlycosidicBondOnTerminalSugarCNPS.add(tmpID);
                                //the first is the overall counter, the second one is specific for this molecule
                                tmpTerminalCircularSugarsWithGlycosidicBondCounter += tmpGlycosidicBondAndTerminalCounter;
                            }
                            if (tmpGlycosidicBondAndNonTerminalCounter > 0) {
                                tmpHasGlycosidicBondOnNonTerminalSugarCounter++;
                                tmpHasGlycosidicBondOnNonTerminalSugarCNPS.add(tmpID);
                            }
                        }
                    }
                    if (tmpHasAnyLinearSugar) {
                        tmpContainsLinearSugarsCounter++;
                        tmpContainsLinearSugarsCNPS.add(tmpID);
                        if (!tmpHasAnyCircularSugar) {
                            tmpHasOnlyLinearSugarsCounter++;
                            tmpHasOnlyLinearSugarsCNPS.add(tmpID);
                        }
                        //terminal and non-terminal
                        List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule);
                        int tmpAllLinearMoietiesCounter;
                        int tmpTerminalCounter;
                        int tmpNonTerminalCounter;
                        tmpAllLinearMoietiesCounter = tmpLinearSugarCandidatesList.size();
                        Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyLinearSugarsMap.get(tmpAllLinearMoietiesCounter);
                        if (Objects.isNull(tmpCurrentValue)) {
                            tmpHowManyMoleculesHaveHowManyLinearSugarsMap.put(tmpAllLinearMoietiesCounter, 1);
                        } else {
                            tmpHowManyMoleculesHaveHowManyLinearSugarsMap.put(tmpAllLinearMoietiesCounter, tmpCurrentValue + 1);
                        }
                        //note:linear moieties that become terminal after removal of a circular moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalLinearMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                        // -1 for the deglycosylated core at the beginning of the list
                        tmpTerminalCounter = tmpRemovedTerminalLinearMoieties.size() - 1 ;
                        tmpNonTerminalCounter = tmpAllLinearMoietiesCounter - tmpTerminalCounter;
                        Assert.assertTrue(tmpNonTerminalCounter >= 0);
                        if (tmpTerminalCounter > 0) {
                            tmpHasTerminalLinearSugarsCounter++;
                            tmpHasTerminalLinearSugarsCNPS.add(tmpID);
                            if (tmpNonTerminalCounter == 0) {
                                tmpHasOnlyTerminalLinearSugarsCounter++;
                                tmpHasOnlyTerminalLinearSugarsCNPS.add(tmpID);
                            }
                        }
                        if (tmpNonTerminalCounter > 0) {
                            tmpHasNonTerminalLinearSugarsCounter++;
                            tmpHasNonTerminalLinearSugarsCNPS.add(tmpID);
                            if (tmpTerminalCounter == 0) {
                                tmpHasOnlyNonTerminalLinearSugarsCounter++;
                                tmpHasOnlyNonTerminalLinearSugarsCNPS.add(tmpID);
                            }
                        }
                    }
                    if (tmpHasAnyCircularSugar && tmpHasAnyLinearSugar) {
                        tmpHasCircularAndLinearSugarsCounter++;
                        tmpHasCircularAndLinearSugarsCNPs.add(tmpID);
                    }
                    IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                    if (tmpDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                        tmpBasicallyASugarCNPS.add(tmpID);
                        int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpMolecule);
                        if (tmpNumberOfMoieties == 1) {
                            tmpBasicallyASingleSugarUnitCounter++;
                            tmpBasicallyASingleSugarUnitCNPs.add(tmpID);
                        }
                    }
                    IAtomContainer tmpCircularDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularSugars(tmpMolecule, true);
                    if (tmpCircularDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyACircularSugarCounter++;
                        tmpBasicallyACircularSugarCNPS.add(tmpID);
                        int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpMolecule);
                        if (tmpNumberOfMoieties == 1) {
                            tmpBasicallyASingleCircularSugarCounter++;
                            tmpBasicallyASingleCircularSugarCNPs.add(tmpID);
                        }
                    }
                    IAtomContainer tmpLinearDeglycosylatedClone = tmpSugarRemovalUtil.removeLinearSugars(tmpMolecule, true);
                    if (tmpLinearDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyALinearSugarCounter++;
                        tmpBasicallyALinearSugarCNPS.add(tmpID);
                        int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfLinearSugars(tmpMolecule);
                        if (tmpNumberOfMoieties == 1) {
                            tmpBasicallyASingleLinearSugarCounter++;
                            tmpBasicallyASingleLinearSugarCNPs.add(tmpID);
                        }
                    }
                }
            } catch (Exception anException) {
                tmpLogger.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                continue;
            }
        }
        //TODO: Outputs in numbers and percentages
        //TODO: Code generic output of one list
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Sugar containing molecules counter: " + tmpSugarContainingMoleculesCounter);
        System.out.println("");
        System.out.println("Circular sugar containing molecules counter: " + tmpContainsCircularSugarsCounter);
        System.out.println("Only circular sugar containing molecules counter: " + tmpHasOnlyCircularSugarsCounter);
        System.out.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalRingsCounter);
        System.out.println("Only terminal circular sugar containing molecules counter: " + tmpHasOnlyTerminalRingsCounter);
        System.out.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalRingsCounter);
        System.out.println("Only non-terminal circular sugar containing molecules counter: " + tmpHasOnlyNonTerminalRingsCounter);
        System.out.println("Circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondCounter);
        System.out.println("Terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnTerminalSugarCounter);
        System.out.println("Non-terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnNonTerminalSugarCounter);
        System.out.println("");
        System.out.println("Detected circular sugar moieties counter: " + tmpCircularSugarsCounter);
        System.out.println("Detected circular sugar moieties that have a glycosidic bond counter: " + tmpCircularSugarsWithGlycosidicBondCounter);
        System.out.println("Detected circular sugar moieties that have a glycosidic bond and are terminal counter: " + tmpTerminalCircularSugarsWithGlycosidicBondCounter);
        System.out.println("");
        System.out.println("Linear sugar containing molecules counter: " + tmpContainsLinearSugarsCounter);
        System.out.println("Only linear sugar containing molecules counter: " + tmpHasOnlyLinearSugarsCounter);
        System.out.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        System.out.println("Only terminal linear sugar containing molecules counter: " + tmpHasOnlyTerminalLinearSugarsCounter);
        System.out.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        System.out.println("Only non-terminal linear sugar containing molecules counter: " + tmpHasOnlyNonTerminalLinearSugarsCounter);
        System.out.println("");
        System.out.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        System.out.println("");
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        System.out.println("Basically a single sugar unit counter: " + tmpBasicallyASingleSugarUnitCounter);
        System.out.println("Basically a circular sugar counter: " + tmpBasicallyACircularSugarCounter);
        System.out.println("Basically a single circular sugar counter: " + tmpBasicallyASingleCircularSugarCounter);
        System.out.println("Basically a linear sugar counter: " + tmpBasicallyALinearSugarCounter);
        System.out.println("Basically a single linear sugar counter: " + tmpBasicallyASingleLinearSugarCounter);
        System.out.println("");
        System.out.println("How many molecules have how many sugars: ");
        for (int i : tmpHowManyMoleculesHaveHowManySugarsMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManySugarsMap.get(i));
        }
        System.out.println("");
        System.out.println("How many molecules have how many circular sugars: ");
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarsMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyCircularSugarsMap.get(i));
        }
        System.out.println("");
        System.out.println("How many molecules have how many linear sugars: ");
        for (int i : tmpHowManyMoleculesHaveHowManyLinearSugarsMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyLinearSugarsMap.get(i));
        }
        tmpCursor.close();
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
        //Note: Here, additional molecules could be added to the list to also test them
        List<String> tmpLinearSugarsList = tmpSugarRemovalUtil.getLinearSugars();
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
     *
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
        System.out.println(tmpSugarRemovalUtil.isGlycosidicBondDetected());
        System.out.println(tmpSugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println(tmpSugarRemovalUtil.getStructuresToKeepMode());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeThreshold());
        System.out.println(tmpSugarRemovalUtil.isNrOfAttachedOxygensIncluded());
        System.out.println(tmpSugarRemovalUtil.getAttachedOxygensToAtomsInRingRatioThreshold());
        System.out.println(tmpSugarRemovalUtil.areLinearSugarsInRingsRemoved());
        System.out.println(tmpSugarRemovalUtil.arePropertiesOfSugarContainingMoleculesSet());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMinSize());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMaxSize());
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
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1OC2=CC(=CC(OC3OC(CO)C(O)C(O)C3O)=C2C4=C1CCC4)C"); //CNP0000001 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //A simple example, the sugar has a glycosidic bond and is not terminal and therefore removed; The resulting
        // disconnected CH3OH is too small to keep and gets cleared away
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
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4"); //CNP0000012 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The sugar ring is not terminal and should not be removed, so the molecule remains unchanged
        Assert.assertEquals("O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
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
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(6);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, one of the phosphate groups is removed because it has only 5 heavy atoms and therefore, the sugar is
        // no longer terminal and also removed
        Assert.assertEquals("O=P(O)(O)OC", tmpSmilesCode);
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(7);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, both phosphate groups are removed because they are too small and nothing remains of the molecule
        Assert.assertEquals("", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        //back to default
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(5);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now, also non-terminal sugars are removed, which leaves two unconnected phosphate groups in this case
        Assert.assertEquals("O=P(O)(O)O.O=P(O)(O)OC", tmpSmilesCode);
        tmpSugarRemovalUtil.setAttachedOxygensToAtomsInRingRatioThreshold(0.7);
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
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //now too, nothing should happen
        Assert.assertEquals("O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CC=C(N)[NH+]=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C%11=CC=[NH+]C=%12NC(NC)CC(C%12%11)CC%10)CO)NCC", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
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
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
        Assert.assertTrue(tmpSugarRemovalUtil.hasLinearSugars(tmpOriginalMolecule));
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Nothing should be removed here because the linear sugar in the macrocycle is detected but not terminal
        Assert.assertEquals("O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
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
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
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
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2OC%10OC(CO)C(O)C(O)C%10O)CO)C(O)C1O"); //CNP0002871 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("OC1CC2CCC3C(CCC4(C)C3CC5OC6(OCC(C)CC6)C(C)C54)C2(C)CC1O", tmpSmilesCode);
        System.out.println();
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMolsAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC4C5=CCC6C7(C)CCC(O)C(C(=O)OC8OC(CO)C(O)C(O)C8O)(C)C7CCC6(C)C5(C)CC3)C(O)C(OC9OC(CO)C(O)C(O)C9O)C2O)C(OC%10OC(CO)C(O)C(O)C%10O)C(O)C1O"); //CNP0005247 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(O)C1(C)C(O)CCC2(C)C1CCC3(C)C2CC=C4C5CC(C)(C)CCC5(C(=O)O)CCC43C", tmpSmilesCode);
        System.out.println();
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMolsAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OC1C(O)C(OC2C3=C(O)C(=CC=C3OC2C(=C)CO)C(=O)C)OC(CO)C1O"); //CNP0032326 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C(C1=CC=C2OC(C(=C)CO)C(O)C2=C1O)C", tmpSmilesCode);
        System.out.println();
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMolsAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"); //CNP0031401 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before)
        Assert.assertEquals("O=C1C=C(OC=2C=C(O)C=C(O)C12)C=3C=CC(O)=C(O)C3", tmpSmilesCode);
        System.out.println();
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMolsAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));

        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2OC3OCC(O)C(OC4OCC(O)C(O)C4O)C3O)OC5C(OC(=O)C67CCC(C)(C)CC7C8=CCC9C%10(C)CC(O)C(OC%11OC(CO)C(O)C(O)C%11O)C(C(=O)O)(C)C%10CCC9(C)C8(CO)CC6)OC(C)C(OC(=O)C=CC%12=CC(OC)=C(OC)C(OC)=C%12)C5OC%13OC(C)C(O)C(O)C%13O)C1O"); //CNP0028122 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //All sugars get removed although some circular sugars only become terminal after the removal of the linear ones
        // (that was a problem before), only one non-terminal sugar remains
        Assert.assertEquals("O=C(OC1C(O)C(O)C(OC(=O)C23CCC(C)(C)CC3C4=CCC5C6(C)CC(O)C(O)C(C(=O)O)(C)C6CCC5(C)C4(CO)CC2)OC1C)C=CC7=CC(OC)=C(OC)C(OC)=C7", tmpSmilesCode);
        System.out.println();
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMolsAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));
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
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
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
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4"); //CNP0032817 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The sugar moiety is NOT connected to the core structure via a glycosidic bond, so it is not removed
        Assert.assertEquals("O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Now that this setting is changed, the sugar moiety is removed
        //note: chemically, the carboxy group should be part of the sugar, not of the core
        Assert.assertEquals("O=COC=1C=C2OC(=CC(=O)C2=C(O)C1O)C=3C=CC(O)=CC3", tmpSmilesCode);

        //another examples for the same thing:
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C([O-])CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"); //CNP0031401 in COCONUTfebruary20
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //The circular sugar moiety is NOT connected to the core structure via a glycosidic bond, so it is not removed
        //also, the removal of linear sugars leaves the circular sugar untouched
        Assert.assertEquals("O=C1C=C(OC=2C1=C(O)C=C(O)C2C3OC(CO)C(O)C(O)C3O)C=4C=CC(O)=C(O)C4", tmpSmilesCode);
        tmpSugarRemovalUtil.setDetectGlycosidicBond(false);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
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
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=CC1(C)C(OC2OC(C(=O)O)C(O)C(OC3OCC(O)C(O)C3O)C2OC4OC(CO)C(O)C(O)C4O)CCC5(C)C6CC=C7C8CC(C)(C)CCC8(C(=O)OC9OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(OC%10OC(CO)C(O)C%10O)C(C)CC)C(C)CC)C(O)C9OC%11OC(C)C(OC%12OCC(O)C(O)C%12O)C(O)C%11O)C(O)CC7(C)C6(C)CCC15"); //CNP0000306 in COCONUTfebruary20
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        //System.out.println(tmpSmilesCode);
        //Only the core structure remains
        Assert.assertEquals("O=CC1(C)C(O)CCC2(C)C1CCC3(C)C2CC=C4C5CC(C)(C)CCC5(C(=O)O)C(O)CC43C", tmpSmilesCode);
        // for illustrative purposes, prints the deglycosylated molecule and all the removed sugar moieties as SMILES strings
        this.printAllMolsAsSmiles(tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpOriginalMolecule, true));
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
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
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
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        //Again, nothing is removed because the linear sugar is non-terminal AND contained in a macrocycle
        Assert.assertEquals("O=C1OC2CC3(OC(CCC3)CC(O)CC)OC(CC4(O)OC(C)(C)CC4C=CCCCCCC(O)(C)C(O)C(OC5OC(C)C(N(C)C)CC5)C(O)C(C)C(O)C(O)(C=C1)C)C2C", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
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
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C"); //CNP0000225 in COCOCNUTfebruary20
        List<IAtomContainer> tmpRemovedSugars = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpOriginalMolecule, true);
        for (IAtomContainer tmpSugar : tmpRemovedSugars) {
            System.out.println(tmpSmiGen.create(tmpSugar));
        }
    }

    /**
     *
     */
    @Test
    public void mariasExampleTest() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator((SmiFlavor.Canonical));
        IAtomContainer tmpOriginalMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpSmilesCode;
        SugarRemovalUtility tmpSugarRemovalUtil = this.getSugarRemovalUtilityV0100DefaultSettings();
        tmpSugarRemovalUtil.setDetectGlycosidicBond(true);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true);
        // a molecule containing 2 circular sugars and 2 linear sugars (sugar acids), in both cases 1 terminal and 1 non-terminal
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C(O)CC(C)(O)CC(=O)OCc1ccc(cc1)OC1C(CO)OC(OC(C(=O)OCc2ccc(OC3OC(CO)CC(O)C3O)cc2)C(O)(CC(C)C)C(=O)OC2CCc3cc4cc(O)c(C)c(O)c4c(O)c3C2=O)C(O)C1O");
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // the terminal linear and circular sugar moieties are removed
        Assert.assertEquals("O=C(OCC1=CC=C(O)C=C1)C(OC2OC(CO)C(OC3=CC=C(C=C3)C)C(O)C2O)C(O)(C(=O)OC4C(=O)C5=C(O)C6=C(O)C(=C(O)C=C6C=C5CC4)C)CC(C)C", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // all 4 sugar moieties are removed, the two non-terminal sugar moieties included
        Assert.assertEquals("O=C1C2=C(O)C3=C(O)C(=C(O)C=C3C=C2CCC1)C.OC1=CC=C(C=C1)C.OC1=CC=C(C=C1)C", tmpSmilesCode);
        // back to default setting
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeLinearSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // only the terminal linear sugar is removed
        Assert.assertEquals("O=C(OCC1=CC=C(OC2OC(CO)CC(O)C2O)C=C1)C(OC3OC(CO)C(OC4=CC=C(C=C4)C)C(O)C3O)C(O)(C(=O)OC5C(=O)C6=C(O)C7=C(O)C(=C(O)C=C7C=C6CC5)C)CC(C)C", tmpSmilesCode);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(false);
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
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(true);
        tmpMoleculeWithoutSugars = tmpSugarRemovalUtil.removeCircularSugars(tmpOriginalMolecule, true);
        tmpSmilesCode = tmpSmiGen.create(tmpMoleculeWithoutSugars);
        System.out.println(tmpSmilesCode);
        // only the terminal circular sugar moiety is removed
        Assert.assertEquals("O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(OC2C(O)C(O)C(OC(C(=O)OCC3=CC=C(O)C=C3)C(O)(C(=O)OC4C(=O)C5=C(O)C6=C(O)C(=C(O)C=C6C=C5CC4)C)CC(C)C)OC2CO)C=C1", tmpSmilesCode);
    }
    //</editor-fold>

    /**
     *
     */
    @Test
    public void classPropertiesTest() throws Exception {
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(true); //default would be false
        tmpSugarRemovalUtil.addLinearSugar("CO");
        List<String> tmpLinearSugarsList = tmpSugarRemovalUtil.getLinearSugars();
        for (String tmpSmiles : tmpLinearSugarsList) {
            System.out.println(tmpSmiles);
        }
        System.out.println("-------------------");
        tmpSugarRemovalUtil.addCircularSugar("O1CCCCCCC1");
        List<String> tmpCircularSugarsList = tmpSugarRemovalUtil.getCircularSugars();
        for (String tmpSmiles : tmpCircularSugarsList) {
            System.out.println(tmpSmiles);
        }
        System.out.println("-------------------");
        System.out.println(tmpSugarRemovalUtil.isGlycosidicBondDetected());
        System.out.println(tmpSugarRemovalUtil.areOnlyTerminalSugarsRemoved());
        System.out.println(tmpSugarRemovalUtil.getStructuresToKeepMode());
        System.out.println(tmpSugarRemovalUtil.getStructureToKeepModeThreshold());
        System.out.println(tmpSugarRemovalUtil.isNrOfAttachedOxygensIncluded());
        System.out.println(tmpSugarRemovalUtil.getAttachedOxygensToAtomsInRingRatioThreshold());
        System.out.println(tmpSugarRemovalUtil.areLinearSugarsInRingsRemoved());
        System.out.println(tmpSugarRemovalUtil.arePropertiesOfSugarContainingMoleculesSet());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMinSize());
        System.out.println(tmpSugarRemovalUtil.getLinearSugarCandidateMaxSize());
        System.out.println(tmpSugarRemovalUtil.areLinearAcidicSugarsDetected());
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
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(true);
        //Interesting case because there is a macrocycle containing a sugar cycle that is not isolated
        tmpOriginalMolecule = tmpSmiPar.parseSmiles("O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO"); //CNP0000509 in COCONUTfebruary20
        List<IAtomContainer> tmpCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpOriginalMolecule);
        for (int i = 0; i < tmpCandidates.size(); i++) {
            IAtomContainer tmpCandidate = tmpCandidates.get(i);
            List<IAtomContainer> tmpToHighlight = new ArrayList<>(1);
            tmpToHighlight.add(tmpCandidate);
            tmpDepictionGenerator.withHighlight(tmpToHighlight, Color.BLUE)
                    .withSize(2000, 2000)
                    .withFillToFit()
                    .depict(tmpOriginalMolecule)
                    .writeTo(tmpOutputFolderPath + File.separator + "Test" + i + ".png");
        }

    }

    //<editor-fold desc="Protected methods">
    /**
     * TODO
     * @return
     */
    protected SugarRemovalUtility getSugarRemovalUtilityV0100DefaultSettings() {
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setDetectGlycosidicBond(false);
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugars(true);
        tmpSugarRemovalUtil.setStructuresToKeepMode(StructuresToKeepMode.HEAVY_ATOM_COUNT);
        tmpSugarRemovalUtil.setStructuresToKeepThreshold(5);
        tmpSugarRemovalUtil.setIncludeNrOfAttachedOxygens(true);
        tmpSugarRemovalUtil.setAttachedOxygensToAtomsInRingRatioThreshold(0.5);
        tmpSugarRemovalUtil.setRemoveLinearSugarsInRing(false);
        tmpSugarRemovalUtil.setLinearSugarCandidateMinSize(4);
        tmpSugarRemovalUtil.setLinearSugarCandidateMaxSize(7);
        tmpSugarRemovalUtil.setDetectLinearAcidicSugars(false);
        tmpSugarRemovalUtil.setPropertyOfSugarContainingMolecules(true);
        return tmpSugarRemovalUtil;
    }
    //</editor-fold>
}
