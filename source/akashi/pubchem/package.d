module akashi.pubchem;

public import akashi.pubchem.assay;
public import akashi.pubchem.compound;
public import akashi.pubchem.conformer3d;
public import akashi.pubchem.protein;
import akashi.pubchem.internal;

import std.algorithm : map;
import std.array : array;
import std.conv : to;
import std.string : join;

Compound[] getProperties(string TYPE)(int[] ids...)
    if (TYPE == "cid" || TYPE == "sid")
{
    return internalGetProperties!TYPE(ids.map!(x => x.to!string).array.join(","));
}

Compound getProperties(string name)
{
    return internalGetProperties!"name"(name)[0];
}

Compound[] getID(string TYPE)(int[] ids...)
    if (TYPE == "cid" || TYPE == "sid")
{
    return internalGetID!TYPE(ids.map!(x => x.to!string).array.join(","));
}

Compound getID(string name)
{
    return internalGetID!"name"(name)[0];
}

Conformer3D[] getConformer3D(string TYPE)(int[] ids...)
    if (TYPE == "cid" || TYPE == "sid")
{
    return internalGetConformer3D!TYPE(ids.map!(x => x.to!string).array.join(","));
}

Conformer3D getConformer3D(string name)
{
    return internalGetConformer3D!"name"(name)[0];
}

string[][] getSynonyms(string TYPE)(int[] ids...)
    if (TYPE == "cid" || TYPE == "sid")
{
    return internalGetSynonyms!TYPE(ids.map!(x => x.to!string).array.join(","));
}

string[] getSynonyms(string name)
{
    return internalGetSynonyms!"name"(name)[0];
}

string[] getDescription(string TYPE)(int[] ids...)
    if (TYPE == "cid" || TYPE == "sid")
{
    return internalGetDescription!TYPE(ids.map!(x => x.to!string).array.join(","));
}

string[] getDescription(string name)
{
    return internalGetDescription!"name"(name);
}

Assay getAssay(int aid)
{
    return internalGetAssay(aid.to!string)[0];
}

AssayResult[] getAssaySummary(string TYPE)(int id)
    if (TYPE == "cid" || TYPE == "sid")
{
    return internalGetAssaySummary!TYPE(id.to!string);
}

Protein getProtein(string accession)
{
    return internalGetProtein(accession)[0];
}

Protein getProteinDetails(string accession)
{
    return internalGetProteinDetails(accession);
}

Assay[] getAssaysByProtein(string accession)
{
    return internalGetAssaysByProtein(accession);
}

Compound[] getCompoundsByAssay(int aid, string cidsType = "all")
{
    return internalGetCompoundsByAssay(aid.to!string, cidsType);
}

Compound[] similaritySearch(int cid, int threshold = 90, int maxRecords = 10)
{
    return internalSimilaritySearch!"cid"(cid.to!string, threshold, maxRecords);
}