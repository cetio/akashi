module akashi.pubchem.internal;

import std.json : JSONValue, parseJSON, JSONType;
import std.conv : to;
import std.uri : encode;
import std.string : assumeUTF, join, lastIndexOf;

import akashi.composer;
import akashi.pubchem.assay;
import akashi.pubchem.conformer3d;
import akashi.pubchem.compound;
import akashi.pubchem.protein;

// Names do not support batch lookup, but CID and SID do.
// TODO: PubMed ID xref support

package:

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

Assay[] internalGetAssay(string str)
{
    Assay[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/assay/aid/"~str~"/summary/JSON"),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "AssaySummaries" !in json || "AssaySummary" !in json["AssaySummaries"])
                throw new Exception("Assay summary is invalid "~json.toString());

            foreach (summary; json["AssaySummaries"]["AssaySummary"].array)
            {
                int aid = intField(summary, "AID");
                if (aid == 0)
                    continue;

                Assay assay = Assay.getOrCreate(aid);
                assay.sourceName = stringField(summary, "SourceName");
                assay.sourceID = stringField(summary, "SourceID");
                assay._name = stringField(summary, "Name");
                assay._description = arrayField(summary, "Description");
                assay._protocol = arrayField(summary, "Protocol");
                assay._comment = arrayField(summary, "Comment");
                assay.numberOfTIDs = intField(summary, "NumberOfTIDs");
                assay.hasScore = boolField(summary, "HasScore");
                assay.method = stringField(summary, "Method");
                assay._targets = assayTargets(summary);
                assay.assayVersion = intField(summary, "Version");
                assay.revision = intField(summary, "Revision");
                assay.lastDataChange = assayDateField(summary, "LastDataChange");
                assay.sidCounts = assayCounts(summary, "SIDCount");
                assay.cidCounts = assayCounts(summary, "CIDCount");
                assay._summaryLoaded = true;
                ret ~= assay;
            }
        },
        null
    );
    return ret;
}

AssayResult[] internalGetAssaySummary(string TYPE)(string str)
{
    AssayResult[] ret;
    string recordNamespace = TYPE == "sid" ? "substance" : "compound";
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/"~recordNamespace~"/"~TYPE~"/"~str~"/assaysummary/JSON"),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "Table" !in json || "Columns" !in json["Table"]
                || "Column" !in json["Table"]["Columns"] || "Row" !in json["Table"])
                return;

            string[] columns = arrayField(json["Table"]["Columns"], "Column");
            foreach (row; json["Table"]["Row"].array)
            {
                if ("Cell" !in row)
                    continue;

                AssayResult result;
                JSONValue[] cells = row["Cell"].array;

                foreach (i, column; columns)
                {
                    string cell = i < cells.length ? jsonText(cells[i]) : null;
                    if (cell.length == 0)
                        continue;

                    if (column == "AID")
                    {
                        int aid = parseInt(cell);
                        if (aid != 0)
                            result.assay = Assay.getOrCreate(aid);
                    }
                    else if (column == "Panel Member ID")
                        result.panelMemberID = parseInt(cell);
                    else if (column == "SID")
                        result.sid = parseInt(cell);
                    else if (column == "CID")
                    {
                        int cid = parseInt(cell);
                        if (cid != 0)
                            result.compound = Compound.getOrCreate(cid);
                    }
                    else if (column == "Activity Outcome")
                        result.outcome = cell;
                    else if (column == "Target Accession")
                        result.targetAccession = cell;
                    else if (column == "Target GeneID")
                        result.targetGeneID = parseInt(cell);
                    else if (column == "Activity Name")
                        result.activityName = cell;
                    else if (column == "Assay Name")
                        result.assayName = cell;
                    else if (column == "Assay Type")
                        result.assayType = cell;
                    else if (column == "PubMed ID")
                        result.pubmedID = parseInt(cell);
                    else if (column == "RNAi")
                        result.rnai = parseBool(cell);
                    else if (column.length >= 14 && column[0 .. 14] == "Activity Value")
                    {
                        result.activityValue = parseDouble(cell);
                        result.activityValueUnit = unitFromColumn(column);
                    }
                    else
                        result.extras[column] = cell;
                }

                if (result.assay is null)
                    continue;

                if (result.compound is null && TYPE == "cid")
                {
                    int cid = parseInt(str);
                    if (cid != 0)
                        result.compound = Compound.getOrCreate(cid);
                }

                if (result.assayName.length > 0 && result.assay._name.length == 0)
                    result.assay._name = result.assayName;

                if (result.targetAccession.length > 0 && !hasTarget(result.assay, result.targetAccession))
                {
                    AssayTarget target;
                    target.accession = result.targetAccession;
                    result.assay._targets ~= target;
                }

                ret ~= result;
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

Assay[] internalGetAssaysByProtein(string accession)
{
    int[] aids = targetAssayIDs("accession", accession);
    if (aids.length == 0)
    {
        Protein protein = internalGetProteinDetails(accession);
        if (protein._geneSymbol.length > 0)
            aids = targetAssayIDs("genesymbol", protein._geneSymbol);
    }

    Assay[] ret;
    foreach (aid; aids)
        ret ~= Assay.getOrCreate(aid);
    return ret;
}

Compound[] internalGetCompoundsByAssay(string str, string cidsType = "all")
{
    Compound[] ret;
    string[string] queryParams;
    if (cidsType.length > 0)
        queryParams["cids_type"] = cidsType;

    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/assay/aid/"~str~"/cids/JSON", queryParams),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "IdentifierList" !in json || "CID" !in json["IdentifierList"])
                return;

            foreach (cidValue; json["IdentifierList"]["CID"].array)
            {
                int cid = valueInt(cidValue);
                if (cid != 0)
                    ret ~= Compound.getOrCreate(cid);
            }
        },
        null
    );
    return ret;
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
            {
                ret ~= Compound.getOrCreate(cid.get!int);
            }
        }, 
        null
    );
    return ret;
}

private:

void parseAtoms(Conformer3D conformer, JSONValue json)
{
    if ("atoms" !in json)
        return;

    JSONValue atoms = json["atoms"];
    int[] aids;
    int[] elements;

    if ("aid" in atoms)
    {
        foreach (aid; atoms["aid"].array)
            aids ~= aid.get!int;
    }

    if ("element" in atoms)
    {
        foreach (elem; atoms["element"].array)
            elements ~= elem.get!int;
    }

    size_t atomCount = aids.length > elements.length ? elements.length : aids.length;
    foreach (size_t i; 0..atomCount)
    {
        Element elem = (elements[i] >= 1 && elements[i] <= 118)
            ? cast(Element)elements[i]
            : Element.H;
        Atom3D atom;
        atom.aid = aids[i];
        atom.element = elem;
        conformer.aidToIndex[aids[i]] = cast(int)i;
        conformer.atoms ~= atom;
    }
}

void parseBonds(Conformer3D conformer, JSONValue json)
{
    if ("bonds" !in json)
        return;

    JSONValue bonds = json["bonds"];
    int[] aid1s;
    int[] aid2s;
    int[] orders;

    if ("aid1" in bonds)
    {
        foreach (aid1; bonds["aid1"].array)
            aid1s ~= aid1.get!int;
    }

    if ("aid2" in bonds)
    {
        foreach (aid2; bonds["aid2"].array)
            aid2s ~= aid2.get!int;
    }

    if ("order" in bonds)
    {
        foreach (order; bonds["order"].array)
            orders ~= order.get!int;
    }

    size_t bondCount = aid1s.length > aid2s.length ? aid2s.length : aid1s.length;
    foreach (i; 0..bondCount)
        conformer.bonds ~= Bond3D(aid1s[i], aid2s[i], i < orders.length ? orders[i] : 1);
}

void parseCoords(Conformer3D conformer, JSONValue json)
{
    if ("coords" !in json || json["coords"].array.length == 0)
        return;

    JSONValue coordSet = json["coords"].array[0];

    if ("conformers" !in coordSet || coordSet["conformers"].array.length == 0)
        return;

    JSONValue conf = coordSet["conformers"].array[0];

    double[] xs, ys, zs;

    if ("x" in conf)
        foreach (x; conf["x"].array)
            xs ~= x.type == JSONType.float_ ? x.get!double : x.get!long.to!double;

    if ("y" in conf)
        foreach (y; conf["y"].array)
            ys ~= y.type == JSONType.float_ ? y.get!double : y.get!long.to!double;

    if ("z" in conf)
        foreach (z; conf["z"].array)
            zs ~= z.type == JSONType.float_ ? z.get!double : z.get!long.to!double;

    if ("aid" in coordSet)
    {
        JSONValue[] coordAids = coordSet["aid"].array;
        foreach (i, aid; coordAids)
        {
            int idx = conformer.indexOf(aid.get!int);
            if (idx >= 0 && i < xs.length && i < ys.length && i < zs.length)
            {
                conformer.atoms[idx].x = xs[i];
                conformer.atoms[idx].y = ys[i];
                conformer.atoms[idx].z = zs[i];
            }
        }
    }

    if ("data" in conf)
    {
        foreach (point; conf["data"].array)
        {
            if ("urn" !in point || "value" !in point)
                continue;

            JSONValue urn = point["urn"];
            JSONValue val = point["value"];

            string label = "label" in urn ? urn["label"].str : null;
            string name = "name" in urn ? urn["name"].str : null;

            if (label == "Conformer" && name == "ID" && "sval" in val)
                conformer.id = val["sval"].str;
            else if (label == "Energy" && "fval" in val)
                conformer.energy = val["fval"].get!double;
            else if (label == "Shape" && name == "Volume" && "fval" in val)
                conformer.volume = val["fval"].get!double;
            else if (label == "Shape" && name == "Self Overlap" && "fval" in val)
                conformer.selfOverlap = val["fval"].get!double;
            else if (label == "Shape" && name == "Multipoles" && "fvec" in val)
            {
                foreach (m; val["fvec"].array)
                {
                    double mval = m.type == JSONType.float_ ? m.get!double : m.get!long.to!double;
                    conformer.multipoles ~= mval;
                }
            }
        }
    }
}

string stringField(JSONValue json, string key)
{
    return key in json ? jsonText(json[key]) : null;
}

string[] arrayField(JSONValue json, string key)
{
    string[] ret;
    if (key !in json)
        return ret;

    foreach (value; json[key].array)
        ret ~= jsonText(value);
    return ret;
}

int intField(JSONValue json, string key)
{
    return key in json ? valueInt(json[key]) : 0;
}

bool boolField(JSONValue json, string key)
{
    if (key !in json)
        return false;

    if (json[key].type == JSONType.true_)
        return true;

    if (json[key].type == JSONType.false_)
        return false;

    return parseBool(jsonText(json[key]));
}

int valueInt(JSONValue value)
{
    switch (value.type)
    {
        case JSONType.integer:
            return value.get!int;
        case JSONType.uinteger:
            return value.get!uint.to!int;
        case JSONType.float_:
            return value.get!double.to!int;
        default:
            return parseInt(jsonText(value));
    }
}

string jsonText(JSONValue value)
{
    switch (value.type)
    {
        case JSONType.string:
            return value.str;
        case JSONType.integer:
            return value.get!long.to!string;
        case JSONType.uinteger:
            return value.get!ulong.to!string;
        case JSONType.float_:
            return value.get!double.to!string;
        case JSONType.true_:
            return "true";
        case JSONType.false_:
            return "false";
        case JSONType.null_:
            return null;
        default:
            return value.toString();
    }
}

int parseInt(string text)
{
    if (text.length == 0)
        return 0;

    try
        return text.to!int;
    catch (Exception)
        return 0;
}

double parseDouble(string text)
{
    if (text.length == 0)
        return double.nan;

    try
        return text.to!double;
    catch (Exception)
        return double.nan;
}

bool parseBool(string text)
{
    return text == "1" || text == "true" || text == "TRUE" || text == "yes" || text == "Yes";
}

string unitFromColumn(string column)
{
    ptrdiff_t start = column.lastIndexOf("[");
    ptrdiff_t stop = column.lastIndexOf("]");
    if (start >= 0 && stop > start + 1)
        return column[start + 1 .. stop];
    return null;
}

AssayCounts assayCounts(JSONValue summary, string prefix)
{
    AssayCounts counts;
    counts.all = intField(summary, prefix~"All");
    counts.active = intField(summary, prefix~"Active");
    counts.inactive = intField(summary, prefix~"Inactive");
    counts.inconclusive = intField(summary, prefix~"Inconclusive");
    counts.unspecified = intField(summary, prefix~"Unspecified");
    counts.probe = intField(summary, prefix~"Probe");
    return counts;
}

AssayDate assayDateField(JSONValue summary, string key)
{
    AssayDate ret;
    if (key !in summary)
        return ret;

    JSONValue value = summary[key];
    ret.year = intField(value, "Year");
    ret.month = intField(value, "Month");
    ret.day = intField(value, "Day");
    return ret;
}

AssayTarget[] assayTargets(JSONValue summary)
{
    AssayTarget[] ret;
    if ("Target" !in summary)
        return ret;

    foreach (targetValue; summary["Target"].array)
    {
        AssayTarget target;
        target.accession = stringField(targetValue, "Accession");
        target.name = stringField(targetValue, "Name");
        ret ~= target;
    }
    return ret;
}

bool hasTarget(Assay assay, string accession)
{
    foreach (target; assay._targets)
    {
        if (target.accession == accession)
            return true;
    }
    return false;
}

int[] targetAssayIDs(string targetType, string value)
{
    int[] ret;
    orchestrator.rateLimit();
    orchestrator.client.get(
        orchestrator.buildURL("/assay/target/"~targetType~"/"~value~"/aids/JSON"),
        (ubyte[] data) {
            JSONValue json = parseJSON(data.assumeUTF);
            if (json.isNull || "IdentifierList" !in json || "AID" !in json["IdentifierList"])
                return;

            foreach (aidValue; json["IdentifierList"]["AID"].array)
            {
                int aid = valueInt(aidValue);
                if (aid != 0)
                    ret ~= aid;
            }
        },
        null
    );
    return ret;
}

JSONValue sectionByHeading(JSONValue json, string heading)
{
    if (json.type == JSONType.object)
    {
        if ("TOCHeading" in json && json["TOCHeading"].str == heading)
            return json;

        if ("Section" in json)
        {
            foreach (section; json["Section"].array)
            {
                JSONValue found = sectionByHeading(section, heading);
                if (!found.isNull)
                    return found;
            }
        }
    }
    else if (json.type == JSONType.array)
    {
        foreach (value; json.array)
        {
            JSONValue found = sectionByHeading(value, heading);
            if (!found.isNull)
                return found;
        }
    }

    return JSONValue.init;
}

string[] sectionStrings(JSONValue record, string heading)
{
    string[] ret;
    JSONValue section = sectionByHeading(record, heading);
    if (section.isNull || "Information" !in section)
        return ret;

    foreach (info; section["Information"].array)
    {
        if ("Value" !in info || "StringWithMarkup" !in info["Value"])
            continue;

        foreach (item; info["Value"]["StringWithMarkup"].array)
        {
            if ("String" in item)
                ret ~= item["String"].str;
        }
    }
    return ret;
}

string sectionString(JSONValue record, string heading)
{
    string[] values = sectionStrings(record, heading);
    return values.length > 0 ? values[0] : null;
}

/*
Similarity

This is a special type of compound namespace input that retrieves CIDs by 2D similarity search. It requires a CID, or a SMILES, InChI, or SDF string in the URL path or POST body (InChI and SDF by POST only). Valid output formats are XML, JSON(P), and ASNT/B.

Example:

https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/2244/cids/XML

Similarity search options are specified via URL arguments:
Option	Type	Meaning	Default
Threshold	integer	minimum Tanimoto score for a hit	90
MaxSeconds	integer	maximum search time in seconds	unlimited
MaxRecords	integer	maximum number of hits	2M
listkey	string	restrict to matches within hits from a prior search	none

Example:

https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/C1=NC2=C(N1)C(=O)N=C(N2)N/cids/XML?Threshold=95&MaxRecords=100
*/