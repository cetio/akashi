module akashi.entrez.eutils;

import akashi.request : RequestClient;

import std.json : JSONValue, parseJSON;
import std.uri : encode;
import std.conv : to;
import std.string : assumeUTF;
import std.array : join;

struct TimeFrame
{
    string from;
    string to;
}

private:

static RequestClient client =
    RequestClient("https://eutils.ncbi.nlm.nih.gov/entrez/eutils", 333);

public:

JSONValue esearch(string DB)(
    string term,
    int retmax = 20,
    int retstart = 0,
    string apiKey = null)
{
    string[string] params = [
        "db": DB,
        "term": term,
        "retmax": retmax.to!string,
        "retstart": retstart.to!string,
        "retmode": "json"
    ];

    if (apiKey !is null && apiKey.length > 0)
        params["api_key"] = apiKey;

    ubyte[] data = client.get("/esearch.fcgi", params);
    if (data is null)
        return JSONValue.init;
    return parseJSON(data.assumeUTF);
}

string efetch(string DB)(string id, string rettype = "abstract", string apiKey = null)
{
    return efetch!DB([id], rettype, apiKey);
}

string efetch(string DB)(string[] ids, string rettype = "abstract", string apiKey = null)
{
    if (ids.length == 0)
        return null;

    ubyte[] data;
    string idParam = ids.join(",");
    if (ids.length > 200 || idParam.length > 2000)
    {
        string postData = "db="~encode(DB)
           ~"&retmode=xml&rettype="~encode(rettype)
           ~"&id="~encode(idParam);
        if (apiKey !is null && apiKey.length > 0)
            postData ~= "&api_key="~encode(apiKey);

        data = client.post("/efetch.fcgi", postData);
    }
    else
    {
        string[string] params = [
            "db": DB,
            "retmode": "xml",
            "rettype": rettype,
            "id": idParam
        ];
        if (apiKey !is null && apiKey.length > 0)
            params["api_key"] = apiKey;

        data = client.get("/efetch.fcgi", params);
    }

    if (data is null)
        return null;
    return cast(string)data.assumeUTF;
}

JSONValue elink(string DBFROM, string DBTO)(
    string id,
    string apiKey = null)
{
    string[string] params = [
        "dbfrom": DBFROM,
        "db": DBTO,
        "id": id,
        "retmode": "json"
    ];
    if (apiKey !is null && apiKey.length > 0)
        params["api_key"] = apiKey;

    ubyte[] data = client.get("/elink.fcgi", params);
    if (data is null)
        return JSONValue.init;
    return parseJSON(data.assumeUTF);
}
