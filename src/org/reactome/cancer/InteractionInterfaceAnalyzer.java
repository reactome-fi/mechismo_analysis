package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;
import org.reactome.funcInt.Protein;
import org.reactome.r3.CosmicAnalyzer;
import org.reactome.r3.HGNCAnalyzer;
import org.reactome.r3.Interactome3dAnalyzer;
import org.reactome.r3.ProteinSequenceHandler;
import org.reactome.r3.ProteinSequenceHandler.Sequence;
import org.reactome.r3.UniProtAnalyzer;
import org.reactome.r3.util.MathUtilities;
import org.reactome.structure.model.InteractionMutationProfile;
import org.reactome.structure.model.InteractionStructure;
import org.reactome.structure.model.PDBUniProtMatch;
import org.reactome.structure.model.ProteinChainInfo;
import org.reactome.structure.model.ProteinMutation;
import org.reactome.structure.model.ProteinMutationProfile;
import org.reactome.structure.model.ResiduleMutationProfile;

/**
 * This class is used to analyze mutation enrichment in the protein interaction interface between
 * two entities.
 * @author wug
 *
 */
public class InteractionInterfaceAnalyzer {
    private static final Logger logger = Logger.getLogger(InteractionInterfaceAnalyzer.class);
    
    // An helper object
    private Interactome3dAnalyzer interactome3dAnalyzer;
    // Cache sequence information for quick performance
    private Map<String, Sequence> accToSeq;
    // For mapping from accession to gene
    private Map<String, String> accessionToGene;
    private Map<String, List<ProteinMutation>> geneToMutations;

    public InteractionInterfaceAnalyzer() {
    }
    
    private void ensurePrecondtions() throws IOException {
        if (interactome3dAnalyzer == null)
            interactome3dAnalyzer = new Interactome3dAnalyzer();
        if (accToSeq == null) {
            logger.info("Loading UniProt accession to sequence...");
            ProteinSequenceHandler sequenceHandler = new ProteinSequenceHandler();
            accToSeq = sequenceHandler.loadSwissProtSequences();
            logger.info("Done.");
        }
        if (accessionToGene == null) {
            logger.info("Loading UniProt accession to gene names...");
            UniProtAnalyzer hgncAnalyzer = new UniProtAnalyzer();
            accessionToGene = hgncAnalyzer.getUniProtAccessionToGeneName();
            logger.info("Done.");
        }
        if (geneToMutations == null || geneToMutations.size() == 0) {
            throw new IllegalStateException("geneToMutations property has not been set!");
        }
    }
    
    public void setGeneToMutations(Map<String, List<ProteinMutation>> geneToMutations) {
        this.geneToMutations = geneToMutations;
    }
    
    public Map<String, List<ProteinMutation>> getGeneToMutations() {
        return this.geneToMutations;
    }
    
    public InteractionMutationProfile analyze(String pdbFileName) throws IOException, StructureException {
        return analyze(new File(pdbFileName));
    }
    
    public InteractionMutationProfile analyze(File pdbFile) throws IOException, StructureException {
        ensurePrecondtions();
        Structure structure = StructureIO.getStructure(pdbFile.getAbsolutePath());
        Map<Chain, List<Integer>> chainToCoordinates = interactome3dAnalyzer.extractContacts(structure);
        String[] tokens = pdbFile.getName().split("-");
        String[] acces = new String[]{tokens[0], tokens[1]};
        Map<Chain, PDBUniProtMatch> chainToMatch = null;
        if (!tokens[2].equals("EXP"))
            chainToMatch = interactome3dAnalyzer.mapCoordinatesToUniProtInPDB(structure,
                    acces,
                    accToSeq,
                    accessionToGene);
        else
            chainToMatch = interactome3dAnalyzer.getMatchForExpStructure(structure,
                    acces,
                    accessionToGene);
        interactome3dAnalyzer.remapCoordinates(chainToMatch, chainToCoordinates);
        InteractionMutationProfile profile = checkInterfaces(chainToCoordinates,
                                                             chainToMatch);
        InteractionStructure intStructure = new InteractionStructure();
        intStructure.setStructure(structure);
        profile.setStructure(intStructure);
        return profile;
    }
    
    private InteractionMutationProfile checkInterfaces(Map<Chain, List<Integer>> chainToContactCoordiantes,
                                                       Map<Chain, PDBUniProtMatch> chainToMatch) {
        // Output from the interface analysis
        // There should be only two entries in these two maps
        if (chainToContactCoordiantes.size() != 2 || chainToMatch.size() != 2)
            throw new IllegalArgumentException("The passed map should contain only two elements!");
        // For two proteins involved in the interaction
        ProteinMutationProfile profile1 = new ProteinMutationProfile();
        ProteinMutationProfile profile2 = new ProteinMutationProfile();
        ProteinMutationProfile current = null;
        CosmicAnalyzer cosmicHelper = new CosmicAnalyzer();
        for (Chain chain : chainToContactCoordiantes.keySet()) {
            // Have to make sure there is a match for the chain
            PDBUniProtMatch match = chainToMatch.get(chain);
            if (match == null)
                throw new IllegalArgumentException("Cannot find a match for chain " + chain);
            if (current == null)
                current = profile1;
            List<ProteinMutation> mutations = geneToMutations.get(match.getGene());
            if (mutations == null) {
                // This should be normal. No need to print out anything here.
//                System.out.println(match.getGene() + " doesn't have an entry in COSMIC!");
                continue;
            }
            // Create a model object to store information
            ProteinChainInfo chainInfo = new ProteinChainInfo();
            List<Integer> contactCoords = chainToContactCoordiantes.get(chain);
            // Use a Protein to store related information
            Protein protein = new Protein();
            protein.setPrimaryAccession(match.getUniprot());
            protein.setShortName(match.getGene());
            String seqKey = "UniProt:" + protein.getPrimaryAccession();
            Sequence sequence = accToSeq.get(seqKey);
            if (sequence == null)
                throw new IllegalStateException("Cannot find sequence for " + protein.getPrimaryAccession());
            protein.setSequence(sequence.getSequence());
            
            chainInfo.setProtein(protein);
            chainInfo.setInterfaceCoordinates(contactCoords);

            current.setChainInfo(chainInfo);
            current.setMutations(geneToMutations.get(protein.getShortName()));

            List<Integer> remappedSeqNumbers = remapChainSeqNumbers(chain, match);
            List<ProteinMutation> mutationsInChain = cosmicHelper.filterEntries(mutations, remappedSeqNumbers);
            if (mutationsInChain.size() == 0) {
                logger.info(match.getGene() + " doesn't have mutations in chain!");
                continue;
            }
            List<ProteinMutation> interfaceMutations = cosmicHelper.filterEntries(mutations, contactCoords);
            if (interfaceMutations.size() == 0) {
                logger.info(match.getGene() + " dosn't have mutations in interface!");
                continue;
            }
            double contactRatio = (double) contactCoords.size() / chain.getAtomGroups().size();
            double pvalue = 1.0;
            if (contactRatio < 1.0d) { // Otherwise, we cannot use binomial test
                pvalue = MathUtilities.calculateBinomialPValue(contactRatio,
                        mutationsInChain.size(),
                        interfaceMutations.size());
                current.setEnrichmentPValue(pvalue);
            }
            ResiduleMutationProfile aaProfile = calculatePValueForPositionMutation(interfaceMutations,
                                                                       mutationsInChain,
                                                                       remappedSeqNumbers);
            current.setMinAAProfile(aaProfile);
        } 
        InteractionMutationProfile intMutProfile = new InteractionMutationProfile();
        intMutProfile.setFirstProteinProfile(profile1);
        intMutProfile.setSecondProteinProfile(profile2);
        return intMutProfile;
    }
    
    /**
     * Calculate p-value for mutation frequency in a single coordinate based on binomial test.
     * @param interfaceMutations
     * @param chainMutations
     * @param chainNumbers
     * @return
     */
    private ResiduleMutationProfile calculatePValueForPositionMutation(List<ProteinMutation> interfaceMutations,
                                                                       List<ProteinMutation> chainMutations,
                                                                       List<Integer> chainNumbers) {
        Map<Integer, Integer> posToMutCount = new HashMap<Integer, Integer>();
        for (ProteinMutation entry : interfaceMutations) {
            Integer count = posToMutCount.get(entry.getCoordinate());
            if (count == null)
                posToMutCount.put(entry.getCoordinate(), 1);
            else
                posToMutCount.put(entry.getCoordinate(), ++count);
        }
        double pvalue = 1.0d; // Largest p-value
        double ratio = 1.0d / chainNumbers.size();
        int minPos = -1;
        for (Integer pos : posToMutCount.keySet()) {
            Integer count = posToMutCount.get(pos);
            double tmpPValue = MathUtilities.calculateBinomialPValue(ratio,
                    chainMutations.size(),
                    count);
            if (tmpPValue < pvalue) {
                pvalue = tmpPValue;
                minPos = pos;
            }
        }
        pvalue = pvalue * posToMutCount.size(); // Perform a Bonferroni correction
        if (pvalue > 1.0d)
            pvalue = 1.0d;
        ResiduleMutationProfile profile = new ResiduleMutationProfile();
        profile.setPvalue(pvalue);
        profile.setCoordinate(minPos);
        return profile;
    }

    
    /**
     * Map coordinates in a chain to the coordinates in UniProt annotation.
     * @param chain
     * @param match
     * @return
     */
    private List<Integer> remapChainSeqNumbers(Chain chain,
                                               PDBUniProtMatch match) {
        List<Integer> list = new ArrayList<Integer>();
        List<Group> groups = chain.getAtomGroups(GroupType.AMINOACID);
        for (Group group : groups) {
            Integer old = group.getResidueNumber().getSeqNum();
            list.add(old + match.getOffset());
        }
        return list;
    }
    
    @Test
    public void testAnalyze() throws Exception {
        Set<String> targetGenes = new HashSet<>();
        
        // Check two genes: EGF and EGFR
        // Cannot find anything with this example
//        targetGenes.add("EGF");
//        targetGenes.add("EGFR");
//        String fileName = "datasets/interactome3d/2017_01/representative/interactions_06/P00533-P01133-EXP-1nql.pdb1-A-0-B-0.pdb";
        
        // Better example
        targetGenes.add("EGFR");
        targetGenes.add("SHC1");
        String fileName = "datasets/interactome3d/2017_01/representative//interactions_06/P00533-P29353-EXP-5czi.pdb1-A-0-B-0.pdb";
        
        CosmicAnalyzer cosmicAnalyzer = new CosmicAnalyzer();
        Map<String, List<ProteinMutation>> geneToCosmicEntries = cosmicAnalyzer.loadMutations(targetGenes);
        System.out.println("Total mutations loaded: " + geneToCosmicEntries.size());
        setGeneToMutations(geneToCosmicEntries);
        InteractionMutationProfile profile = analyze(fileName);
        ProteinMutationProfile proteinProfile = profile.getFirstProteinProfile();
        System.out.println(proteinProfile.getGeneName() + ": " + proteinProfile.getEnrichmentPValue() + ", " + proteinProfile.getMinAAProfile());
        proteinProfile = profile.getSecondProteinProfile();
        System.out.println(proteinProfile.getGeneName() + ": " + proteinProfile.getEnrichmentPValue() + ", " + proteinProfile.getMinAAProfile());
    }

}
