module akashi.pubchem.internal.parse;

import std.json : JSONValue, JSONType;
import std.conv : to;
import std.string : lastIndexOf;

import akashi.pubchem.conformer3d;

package:

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
    foreach (size_t i; 0 .. atomCount)
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
    foreach (i; 0 .. bondCount)
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

    double[] xs;
    double[] ys;
    double[] zs;

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
    JSONValue section = sectionByHeading(record, heading);
    return sectionValueStrings(section);
}

string sectionString(JSONValue record, string heading)
{
    string[] values = sectionStrings(record, heading);
    return values.length > 0 ? values[0] : null;
}

string[][string] sectionChildStrings(JSONValue record, string heading)
{
    string[][string] ret;
    JSONValue section = sectionByHeading(record, heading);
    if (section.isNull || "Section" !in section)
        return ret;

    foreach (child; section["Section"].array)
    {
        if ("TOCHeading" !in child)
            continue;

        string[] values = sectionValueStrings(child);
        if (values.length > 0)
            ret[child["TOCHeading"].str] = values;
    }
    return ret;
}

string[] sectionNamedStrings(JSONValue record, string heading, string name)
{
    JSONValue section = sectionByHeading(record, heading);
    if (section.isNull || "Information" !in section)
        return null;

    foreach (info; section["Information"].array)
    {
        if ("Name" !in info || info["Name"].str != name || "Value" !in info)
            continue;

        return valueStrings(info["Value"]);
    }

    return null;
}

private:

string[] sectionValueStrings(JSONValue section)
{
    string[] ret;
    if (section.isNull || "Information" !in section)
        return ret;

    foreach (info; section["Information"].array)
    {
        if ("Value" !in info)
            continue;

        ret ~= valueStrings(info["Value"]);
    }
    return ret;
}

string[] valueStrings(JSONValue value)
{
    string[] ret;
    if ("StringWithMarkup" !in value)
        return ret;

    foreach (item; value["StringWithMarkup"].array)
    {
        if ("String" in item)
            ret ~= item["String"].str;
    }
    return ret;
}
