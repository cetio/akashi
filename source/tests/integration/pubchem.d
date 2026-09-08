module tests.integration.pubchem;

version(pubchem)
{
    import akashi.pubchem : Compound, Conformer3D, Gene, Protein, getGene, getGenesBySynonym, getID,
        getProperties, getProtein, getProteinDetails, similaritySearch;

    import core.thread : Thread;
    import core.time : dur;
    import std.algorithm : canFind;
    import std.string : indexOf, toLower;

private:
    Compound compound;
    Compound identifiedCompound;
    Compound[] similarCompounds;
    Conformer3D conformer;
    Protein protein;
    Gene gene;
    Gene[] synonymGenes;
    string[] synonyms;
    string description;

    static this()
    {
        compound = getProperties("caffeine");
        synonyms = compound.synonyms;
        description = compound.description;
        conformer = compound.conformer3D;

        Thread.sleep(dur!"seconds"(1));

        identifiedCompound = getID("caffeine");
        similarCompounds = similaritySearch(2519, 99, 5);
        protein = getProtein("P23219");
        protein = getProteinDetails("P23219");

        Thread.sleep(dur!"seconds"(1));

        gene = getGene(7157);
        synonymGenes = getGenesBySynonym("P53");
        assert(gene.description.length > 0);
    }

    unittest
    {
        assert(compound.cid == 2519);
        assert(compound.name == "Caffeine");
    }

    unittest
    {
        assert(compound.properties.formula == "C8H10N4O2");
        assert(compound.properties.weight > 190 && compound.properties.weight < 200);
    }

    unittest
    {
        assert(identifiedCompound.cid == 2519);
        assert(identifiedCompound.sids.length > 0);
    }

    unittest
    {
        assert(synonyms.length > 0);
        assert(synonyms.canFind!(synonym => synonym.toLower == "caffeine"));
    }

    unittest
    {
        assert(description.length > 50);
    }

    unittest
    {
        assert(similarCompounds.length > 0 && similarCompounds.length <= 5);
        assert(similarCompounds.canFind!(candidate => candidate.cid == 2519));
    }

    unittest
    {
        assert(conformer.cid == 2519);
        assert(conformer.isValid);
        assert(conformer.atoms.length > 10);
        assert(conformer.bonds.length > 10);
        assert(conformer.indexOf(conformer.atoms[0].aid) == 0);
    }

    unittest
    {
        assert(protein.accession == "P23219");
        assert(protein.name.length > 0);
        assert(protein.description.length > 0);
        assert(protein.externalURL.indexOf("ncbi.nlm.nih.gov") >= 0);
    }

    unittest
    {
        assert(gene.geneID == 7157);
        assert(gene.symbol == "TP53");
        assert(gene.name.length > 0);
        assert(gene.taxonomy.length > 0);
    }

    unittest
    {
        assert(synonymGenes.length > 0);
        assert(gene.identifiers.length > 0);
    }
}
