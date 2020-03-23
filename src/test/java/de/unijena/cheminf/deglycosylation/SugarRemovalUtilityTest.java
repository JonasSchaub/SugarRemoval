/**
 * TODO: Add license header
 */
package de.unijena.cheminf.deglycosylation;

/**
 * TODO:
 * - Add more examples
 * - play around with settings
 */

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

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
}
