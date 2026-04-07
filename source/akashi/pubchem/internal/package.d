module akashi.pubchem.internal;

import std.json : JSONValue, parseJSON;
import std.conv : to;
import std.uri : encode;
import std.string : assumeUTF, join, replace;

import akashi.composer;
import akashi.pubchem.bio.gene;
import akashi.pubchem.bio.protein;
import akashi.pubchem.conformer3d;
import akashi.pubchem.compound;
import akashi.pubchem.internal.parse;

// Names do not support batch lookup, but CID and SID do.
// TODO: PubMed ID xref support

package(akashi.pubchem):

static Orchestrator orchestrator = Orchestrator("https://pubchem.ncbi.nlm.nih.gov/rest/pug", 200);
static Orchestrator viewOrchestrator = Orchestrator("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view", 200);

/// Fetches `Compound[]` with `cid` and `properties` populated by `/compound/%{TYPE}/%{ids}/property/`
Compound[] internalGetProperties(string TYPE)(string str)
{
    string[] props = [
        "Title", "SMILES", "IUPACName", "InChI",
        "MolecularFormula", "MolecularWeight", "ExactMass", "Charge", "TPSA", "XLogP"
    ];

    Compound[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/compound/"~TYPE~"/"~str~"/property/"~props.join(",")~"/JSON"), 
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "PropertyTable" !in json || "Properties" !in json["PropertyTable"])
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
        }, 
        null
    );
    return ret;
}

/// Fetches `Compound` with `cid` and `sids` populated by `/compound/%{TYPE}/%{str}/sids/JSON`
Compound[] internalGetID(string TYPE)(string str)
{
    Compound[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/compound/"~TYPE~"/"~str~"/sids/JSON"), 
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "InformationList" !in json || "Information" !in json["InformationList"])
                throw new Exception("ID list is invalid "~json.toString());

            foreach (ids; json["InformationList"]["Information"].array)
            {
                Compound compound = Compound.getOrCreate(ids["CID"].get!int);
                foreach (sid; ids["SID"].array)
                    compound.sids ~= sid.get!int;
                ret ~= compound;
            }
        }, 
        null
    );
    return ret;
}

Conformer3D[] internalGetConformer3D(string TYPE)(string str)
{
    Conformer3D[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/compound/"~TYPE~"/"~str~"/JSON", ["record_type": "3d"]), 
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "PC_Compounds" !in json || json["PC_Compounds"].array.length == 0)
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
        }, 
        null
    );
    return ret;
}

string[][] internalGetSynonyms(string TYPE)(string str)
{
    string[][] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/compound/"~TYPE~"/"~str~"/synonyms/JSON"), 
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "InformationList" !in json || "Information" !in json["InformationList"])
                return;

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
        }, 
        null
    );
    return ret;
}

string[] internalGetDescription(string TYPE)(string str)
{
    string[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/compound/"~TYPE~"/"~str~"/description/JSON"), 
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "InformationList" !in json || "Information" !in json["InformationList"])
                return;

            foreach (info; json["InformationList"]["Information"].array)
            {
                if ("Description" in info)
                    ret ~= info["Description"].str;
            }
        }, 
        null
    );
    return ret;
}

Protein[] internalGetProtein(string str)
{
    Protein[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/protein/accession/"~str~"/summary/JSON"),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "ProteinSummaries" !in json || "ProteinSummary" !in json["ProteinSummaries"])
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
        },
        null
    );
    return ret;
}

Protein internalGetProteinDetails(string accession)
{
    Protein protein = Protein.getOrCreate(accession);
    viewOrchestrator.rateLimit();
    viewOrchestrator.client.get(
        viewOrchestrator.buildURL("/data/protein/"~accession~"/JSON"),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "Record" !in json)
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
        },
        null
    );
    return protein;
}

Gene[] internalGetGeneSummary(string path)
{
    Gene[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL(path),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "GeneSummaries" !in json || "GeneSummary" !in json["GeneSummaries"])
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
        },
        null
    );
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
    viewOrchestrator.rateLimit();
    viewOrchestrator.client.get(
        viewOrchestrator.buildURL("/data/gene/"~geneID.to!string~"/JSON"),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "Record" !in json)
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

            foreach (i; 0 .. isoformCount)
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
        },
        null
    );
    return gene;
}

Compound[] internalSimilaritySearch(string TYPE)(string str, int threshold = 90, int maxRecords = 2_000_000)
{
    Compound[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL(
            "/compound/fastsimilarity_2d/"~TYPE~"/"~str~"/cids/JSON",
            ["Threshold": threshold.to!string, "MaxRecords": maxRecords.to!string]
        ),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "IdentifierList" !in json || "CID" !in json["IdentifierList"])
                throw new Exception("Similarity search results are invalid "~json.toString());

            foreach (cid; json["IdentifierList"]["CID"].array)
                ret ~= Compound.getOrCreate(cid.get!int);
        },
        null
    );
    return ret;
}