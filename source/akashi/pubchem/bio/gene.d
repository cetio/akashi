module akashi.pubchem.bio.gene;

import akashi.pubchem.bio : getGene, getGeneDetails, getGenesBySynonym;
import akashi.pubchem.bio.protein;

struct GeneIsoform
{
    string name;
    string uniprotID;
    string refSeqAccession;
}

class Gene
{
    static Gene[int] registry;

    static Gene getOrCreate(int geneID)
    {
        assert(geneID != 0, "GeneID must not be zero");
        if (auto p = geneID in registry)
            return *p;

        Gene gene = new Gene(geneID);
        registry[geneID] = gene;
        return gene;
    }

package(akashi.pubchem):
    this(int geneID)
    {
        assert(geneID != 0, "GeneID must not be zero");
        this.geneID = geneID;
    }

    string _symbol;
    string _name;
    string _taxonomy;
    string[] _synonyms;
    string _description;
    string[][string] _identifiers;
    string[] _orthologs;
    string[] _proteinFunctions;
    string[] _proteinAccessions;
    GeneIsoform[] _proteinIsoforms;
    string[] _diseases;
    string[] _pathways;
    bool _summaryLoaded;
    bool _detailsLoaded;

public:
    int geneID;
    int taxonomyID;
    string externalURL;

    static Gene byID(int geneID)
    {
        return getGene(geneID);
    }

    static Gene bySymbol(string symbol)
    {
        return getGene(symbol);
    }

    static Gene bySymbol(string symbol, int taxonomyID)
    {
        return getGene(symbol, taxonomyID);
    }

    static Gene bySymbol(string symbol, string taxonomy)
    {
        return getGene(symbol, taxonomy);
    }

    static Gene[] bySynonym(string synonym)
    {
        return getGenesBySynonym(synonym);
    }

    ref string symbol()
    {
        ensureSummary();
        return _symbol;
    }

    ref string name()
    {
        ensureSummary();
        return _name;
    }

    ref string taxonomy()
    {
        ensureSummary();
        return _taxonomy;
    }

    string[] synonyms()
    {
        ensureSummary();
        return _synonyms;
    }

    ref string description()
    {
        ensureDetails();
        return _description;
    }

    string[][string] identifiers()
    {
        ensureDetails();
        return _identifiers;
    }

    string[] orthologs()
    {
        ensureDetails();
        return _orthologs;
    }

    string[] proteinFunctions()
    {
        ensureDetails();
        return _proteinFunctions;
    }

    string[] proteinAccessions()
    {
        ensureDetails();
        return _proteinAccessions;
    }

    GeneIsoform[] proteinIsoforms()
    {
        ensureDetails();
        return _proteinIsoforms;
    }

    Protein[] proteins()
    {
        Protein[] ret;
        bool[string] seen;

        foreach (accession; proteinAccessions())
        {
            if (accession.length == 0 || accession in seen)
                continue;

            seen[accession] = true;
            ret ~= Protein.getOrCreate(accession);
        }
        return ret;
    }

    string[] diseases()
    {
        ensureDetails();
        return _diseases;
    }

    string[] pathways()
    {
        ensureDetails();
        return _pathways;
    }

private:
    void ensureSummary()
    {
        if (!_summaryLoaded)
            getGene(geneID);
    }

    void ensureDetails()
    {
        ensureSummary();
        if (!_detailsLoaded)
            getGeneDetails(geneID);
    }
}
