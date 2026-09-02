module akashi.pubchem.internal;

import akashi.pubchem.bio.gene;
import akashi.pubchem.bio.protein;
import akashi.pubchem.conformer3d;
import akashi.pubchem.compound;
import akashi.pubchem.internal.parse;
import akashi.request : RequestClient;

import std.json : JSONValue, parseJSON;
import std.conv : to;
import std.uri : encode;
import std.string : assumeUTF, join, replace;

// Names do not support batch lookup, but CID and SID do.
// TODO: PubMed ID xref support

private:

JSONValue fetchJSON(ref RequestClient client, string path, string[string] query = null)
{
    ubyte[] data = client.get(path, query);
    if (data is null)
        return JSONValue.init;
    return parseJSON(data.assumeUTF);
}

package(akashi.pubchem):

static RequestClient client = RequestClient("https://pubchem.ncbi.nlm.nih.gov/rest/pug", 200);
static RequestClient viewClient = RequestClient("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view", 200);

/// Fetches `Compound[]` with `cid` and `properties` populated by `/compound/%{TYPE}/%{ids}/property/`
Compound[] internalGetProperties(string TYPE)(string str)
{
    string[] propertyNames = [
        "Title",
        "SMILES",
        "IUPACName",
        "InChI",
        "MolecularFormula",
        "MolecularWeight",
        "ExactMass",
        "Charge",
        "TPSA",
        "XLogP"
    ];

    Compound[] ret;
    JSONValue json = fetchJSON(
        client,
        "/compound/"~TYPE~"/"~str~"/property/"~propertyNames.join(",")~"/JSON"
    );
    if (json.isNull)
        return ret;

    if ("PropertyTable" !in json || "Properties" !in json["PropertyTable"])
        throw new Exception("Properties are invalid "~json.toString());

    foreach (props; json["PropertyTable"]["Properties"].array)
    {
        Compound compound = Compound.getOrCreate(props["CID"].get!int);
        compound.properties = Properties(
            "Title" in props ? props["Title"].str : null,
            "SMILES" in props ? props["SMILES"].str : null,
            "IUPACName" in props ? props["IUPACName"].str : null,
            "InChI" in props ? props["InChI"].str : null,

            "MolecularFormula" in props ? props["MolecularFormula"].str : null,
            "MolecularWeight" in props ? props["MolecularWeight"].str.to!double : double.nan,
            "ExactMass" in props ? props["ExactMass"].str.to!double : double.nan,
            "Charge" in props ? props["Charge"].get!int : 0,
            "TPSA" in props ? props["TPSA"].get!double : double.nan,
            "XLogP" in props ? props["XLogP"].get!double : double.nan
        );
        ret ~= compound;
    }
    return ret;
}

/// Fetches `Compound` with `cid` and `sids` populated by `/compound/%{TYPE}/%{str}/sids/JSON`
Compound[] internalGetID(string TYPE)(string str)
{
    Compound[] ret;
    JSONValue json = fetchJSON(client, "/compound/"~TYPE~"/"~str~"/sids/JSON");
    if (json.isNull)
        return ret;

    if ("InformationList" !in json || "Information" !in json["InformationList"])
        throw new Exception("ID list is invalid "~json.toString());

    foreach (ids; json["InformationList"]["Information"].array)
    {
        Compound compound = Compound.getOrCreate(ids["CID"].get!int);
        foreach (sid; ids["SID"].array)
            compound.sids ~= sid.get!int;
        ret ~= compound;
    }
    return ret;
}

Conformer3D[] internalGetConformer3D(string TYPE)(string str)
{
    Conformer3D[] ret;
    JSONValue json = fetchJSON(
        client,
        "/compound/"~TYPE~"/"~str~"/JSON",
        ["record_type": "3d"]
    );
    if (json.isNull)
        return ret;

    if ("PC_Compounds" !in json || json["PC_Compounds"].array.length == 0)
        throw new Exception("Conformer 3D not found for "~str);

    foreach (pc; json["PC_Compounds"].array)
    {
        Conformer3D conformer = new Conformer3D();
        conformer.cid = pc["id"]["id"]["cid"].get!int;
        conformer.parseAtoms(pc);
        conformer.parseBonds(pc);
        conformer.parseCoords(pc);
        ret ~= conformer;
    }
    return ret;
}

string[][] internalGetSynonyms(string TYPE)(string str)
{
    string[][] ret;
    JSONValue json = fetchJSON(client, "/compound/"~TYPE~"/"~str~"/synonyms/JSON");
    if (json.isNull || "InformationList" !in json || "Information" !in json["InformationList"])
        return ret;

    foreach (info; json["InformationList"]["Information"].array)
    {
        string[] synonyms;
        if ("Synonym" in info)
        {
            foreach (syn; info["Synonym"].array)
                synonyms ~= syn.str;
        }
        ret ~= synonyms;
    }
    return ret;
}

string[] internalGetDescription(string TYPE)(string str)
{
    string[] ret;
    JSONValue json = fetchJSON(client, "/compound/"~TYPE~"/"~str~"/description/JSON");
    if (json.isNull || "InformationList" !in json || "Information" !in json["InformationList"])
        return ret;

    foreach (info; json["InformationList"]["Information"].array)
    {
        if ("Description" in info)
            ret ~= info["Description"].str;
    }
    return ret;
}

Protein[] internalGetProtein(string str)
{
    Protein[] ret;
    JSONValue json = fetchJSON(client, "/protein/accession/"~str~"/summary/JSON");
    if (json.isNull)
        return ret;

    if ("ProteinSummaries" !in json || "ProteinSummary" !in json["ProteinSummaries"])
        throw new Exception("Protein summary is invalid "~json.toString());

    foreach (summary; json["ProteinSummaries"]["ProteinSummary"].array)
    {
        string accession = stringField(summary, "ProteinAccession");
        if (accession.length == 0)
            continue;

        Protein protein = Protein.getOrCreate(accession);
        protein._name = stringField(summary, "Name");
        protein.taxonomyID = intField(summary, "TaxonomyID");
        protein._taxonomy = stringField(summary, "Taxonomy");
        protein._synonyms = arrayField(summary, "Synonym");
        protein.externalURL = "https://www.ncbi.nlm.nih.gov/protein/"~accession;
        protein._summaryLoaded = true;
        ret ~= protein;
    }
    return ret;
}

Protein internalGetProteinDetails(string accession)
{
    Protein protein = Protein.getOrCreate(accession);
    JSONValue json = fetchJSON(viewClient, "/data/protein/"~accession~"/JSON");
    if (json.isNull)
        return protein;

    if ("Record" !in json)
        throw new Exception("Protein details are invalid "~json.toString());

    JSONValue record = json["Record"];
    if (protein._name.length == 0)
        protein._name = stringField(record, "RecordTitle");

    if ("RecordExternalURL" in record)
        protein.externalURL = record["RecordExternalURL"].str;

    string taxonomy = sectionString(record, "Taxonomy");
    if (taxonomy.length > 0)
        protein._taxonomy = taxonomy;

    string[] synonyms = sectionStrings(record, "Synonyms");
    if (synonyms.length > 0)
        protein._synonyms = synonyms;

    protein._description = sectionString(record, "Record Description");
    if (protein._description is null)
        protein._description = "";

    protein._geneSymbol = sectionString(record, "Encoding Gene");
    if (protein._geneSymbol is null)
        protein._geneSymbol = "";

    protein._refSeqAccessions = sectionStrings(record, "RefSeq Accession");
    protein._summaryLoaded = true;
    protein._detailsLoaded = true;
    return protein;
}

Gene[] internalGetGeneSummary(string path)
{
    Gene[] ret;
    JSONValue json = fetchJSON(client, path);
    if (json.isNull)
        return ret;

    if ("GeneSummaries" !in json || "GeneSummary" !in json["GeneSummaries"])
        throw new Exception("Gene summary is invalid "~json.toString());

    foreach (summary; json["GeneSummaries"]["GeneSummary"].array)
    {
        int geneID = intField(summary, "GeneID");
        if (geneID == 0)
            continue;

        Gene gene = Gene.getOrCreate(geneID);
        gene._symbol = stringField(summary, "Symbol");
        gene._name = stringField(summary, "Name");
        gene.taxonomyID = intField(summary, "TaxonomyID");
        gene._taxonomy = stringField(summary, "Taxonomy");
        gene._description = stringField(summary, "Description");
        gene._synonyms = arrayField(summary, "Synonym");
        gene.externalURL = "https://www.ncbi.nlm.nih.gov/gene/"~gene.geneID.to!string;
        gene._summaryLoaded = true;
        ret ~= gene;
    }
    return ret;
}

Gene[] internalGetGeneByID(string str)
{
    return internalGetGeneSummary("/gene/geneid/"~str~"/summary/JSON");
}

Gene[] internalGetGeneBySymbol(string symbol, string taxonomy = null)
{
    string path = "/gene/genesymbol/"~encode(symbol).replace("+", "%20");
    if (taxonomy.length > 0)
        path ~= "/"~encode(taxonomy).replace("+", "%20");
    return internalGetGeneSummary(path~"/summary/JSON");
}

Gene[] internalGetGeneBySynonym(string synonym)
{
    return internalGetGeneSummary("/gene/synonym/"~encode(synonym).replace("+", "%20")~"/summary/JSON");
}

Gene internalGetGeneDetails(int geneID)
{
    Gene gene = Gene.getOrCreate(geneID);
    JSONValue json = fetchJSON(viewClient, "/data/gene/"~geneID.to!string~"/JSON");
    if (json.isNull)
        return gene;

    if ("Record" !in json)
        throw new Exception("Gene details are invalid "~json.toString());

    JSONValue record = json["Record"];
    if ("RecordExternalURL" in record)
        gene.externalURL = record["RecordExternalURL"].str;

    if (gene._name.length == 0)
        gene._name = stringField(record, "RecordTitle");

    string symbol = sectionString(record, "Symbol");
    if (symbol.length > 0)
        gene._symbol = symbol;

    string taxonomy = sectionString(record, "Taxonomy");
    if (taxonomy.length > 0)
        gene._taxonomy = taxonomy;

    string[] synonyms = sectionStrings(record, "Synonyms");
    if (synonyms.length > 0)
        gene._synonyms = synonyms;

    string description = sectionString(record, "Record Description");
    if (description !is null)
        gene._description = description;
    else if (gene._description is null)
        gene._description = "";

    gene._identifiers = sectionChildStrings(record, "Other Identifiers");
    gene._orthologs = sectionStrings(record, "Orthologous Genes");
    gene._proteinFunctions = sectionStrings(record, "Protein Function");

    string[] proteinAccessions = sectionStrings(record, "Protein Targets");
    if (proteinAccessions.length > 0)
        gene._proteinAccessions = proteinAccessions;
    else if (gene._proteinAccessions.length == 0)
        gene._proteinAccessions = sectionNamedStrings(record, "Protein Isoforms", "RefSeq Accession");

    gene._proteinIsoforms.length = 0;
    string[] isoformNames = sectionNamedStrings(record, "Protein Isoforms", "Isoform");
    string[] isoformUniprotIDs = sectionNamedStrings(record, "Protein Isoforms", "UniProt ID");
    string[] isoformRefSeqAccessions = sectionNamedStrings(record, "Protein Isoforms", "RefSeq Accession");
    size_t isoformCount = isoformNames.length;
    if (isoformUniprotIDs.length > isoformCount)
        isoformCount = isoformUniprotIDs.length;

    if (isoformRefSeqAccessions.length > isoformCount)
        isoformCount = isoformRefSeqAccessions.length;

    foreach (i; 0..isoformCount)
    {
        GeneIsoform isoform;
        if (i < isoformNames.length)
            isoform.name = isoformNames[i];

        if (i < isoformUniprotIDs.length)
            isoform.uniprotID = isoformUniprotIDs[i];

        if (i < isoformRefSeqAccessions.length)
            isoform.refSeqAccession = isoformRefSeqAccessions[i];

        if (isoform.name.length > 0 || isoform.uniprotID.length > 0 || isoform.refSeqAccession.length > 0)
            gene._proteinIsoforms ~= isoform;
    }

    string[] diseases;
    diseases ~= sectionStrings(record, "GHR Health Conditions");
    diseases ~= sectionStrings(record, "OMIM Phenotypes");
    diseases ~= sectionStrings(record, "MedGen Diseases");
    gene._diseases = diseases;
    gene._pathways = sectionStrings(record, "Pathways");
    gene._summaryLoaded = true;
    gene._detailsLoaded = true;
    return gene;
}

Compound[] internalSimilaritySearch(string TYPE)(string str, int threshold = 90, int maxRecords = 2_000_000)
{
    Compound[] ret;
    JSONValue json = fetchJSON(
        client,
        "/compound/fastsimilarity_2d/"~TYPE~"/"~str~"/cids/JSON",
        ["Threshold": threshold.to!string, "MaxRecords": maxRecords.to!string]
    );
    if (json.isNull)
        return ret;

    if ("IdentifierList" !in json || "CID" !in json["IdentifierList"])
        throw new Exception("Similarity search results are invalid "~json.toString());

    foreach (cid; json["IdentifierList"]["CID"].array)
        ret ~= Compound.getOrCreate(cid.get!int);
    return ret;
}