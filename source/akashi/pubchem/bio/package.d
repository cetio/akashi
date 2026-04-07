module akashi.pubchem.bio;

public import akashi.pubchem.bio.gene;
public import akashi.pubchem.bio.protein;

import akashi.pubchem.internal;

import std.conv : to;

Protein getProtein(string accession)
{
    return internalGetProtein(accession)[0];
}

Protein getProteinDetails(string accession)
{
    getProtein(accession);
    return internalGetProteinDetails(accession);
}

Gene getGene(int geneID)
{
    return internalGetGeneByID(geneID.to!string)[0];
}

Gene getGene(string symbol)
{
    return internalGetGeneBySymbol(symbol)[0];
}

Gene getGene(string symbol, int taxonomyID)
{
    return internalGetGeneBySymbol(symbol, taxonomyID.to!string)[0];
}

Gene getGene(string symbol, string taxonomy)
{
    return internalGetGeneBySymbol(symbol, taxonomy)[0];
}

Gene[] getGenesBySynonym(string synonym)
{
    return internalGetGeneBySynonym(synonym);
}

Gene getGeneDetails(int geneID)
{
    getGene(geneID);
    return internalGetGeneDetails(geneID);
}
