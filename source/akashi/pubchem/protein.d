module akashi.pubchem.protein;

import akashi.pubchem : getProtein, getProteinDetails, getAssaysByProtein;
import akashi.pubchem.assay;

class Protein
{
    static Protein[string] registry;

    static Protein getOrCreate(string accession)
    {
        assert(accession.length > 0, "Protein accession must not be empty");
        if (auto p = accession in registry)
            return *p;

        Protein protein = new Protein(accession);
        registry[accession] = protein;
        return protein;
    }

package:
    this(string accession)
    {
        assert(accession.length > 0, "Protein accession must not be empty");
        this.accession = accession;
    }

    string _name;
    string _taxonomy;
    string[] _synonyms;
    string _description;
    string _geneSymbol;
    string[] _refSeqAccessions;
    bool _summaryLoaded;
    bool _detailsLoaded;
    Assay[] _assays;
    bool _assaysLoaded;

public:
    string accession;
    int taxonomyID;
    string externalURL;

    static Protein byAccession(string accession)
    {
        return getProtein(accession);
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

    ref string geneSymbol()
    {
        ensureDetails();
        return _geneSymbol;
    }

    string[] refSeqAccessions()
    {
        ensureDetails();
        return _refSeqAccessions;
    }

    Assay[] assays()
    {
        if (!_assaysLoaded)
        {
            _assays = getAssaysByProtein(accession);
            _assaysLoaded = true;
        }

        return _assays;
    }

private:
    void ensureSummary()
    {
        if (!_summaryLoaded)
            getProtein(accession);
    }

    void ensureDetails()
    {
        ensureSummary();
        if (!_detailsLoaded)
            getProteinDetails(accession);
    }
}
