module akashi.page;

import std.array : join;
import std.algorithm : canFind, map;

import akashi.text.ast;
import akashi.text.wiki : parseWikitext;
import akashi.text.xml : parseXml;
import akashi.pubchem.compound : Compound;

struct AkashiConfig
{
    string[] enabledSources = ["wikipedia", "psychonaut"];
}

__gshared AkashiConfig config;

class Page
{
package:
    string _raw;
    Document _doc;
    bool _parsed;
    void delegate(Page) _fetchContent;

    this() { }

    this(
        string title,
        string source,
        string url,
        void delegate(Page) fetchContent
    )
    {
        this.title = title;
        this.source = source;
        this.url = url;
        this._fetchContent = fetchContent;
    }

public:
    string title;
    string source;
    string url;

    static Page fromRaw(string title, string source, string rawContent)
    {
        Page p = new Page();
        p.title = title;
        p.source = source;
        p._raw = rawContent;
        return p;
    }

    /// Access the raw source text, fetching lazily if needed.
    ref string raw()
    {
        if (_raw is null && _fetchContent !is null)
        {
            _fetchContent(this);
            _fetchContent = null;
        }
        return _raw;
    }

    /// The parsed AST Document. Lazily parsed on first access.
    /// Non-renderable node types (Templates, References, Comments,
    /// Categories, Images) are dropped after parsing.
    ref Document document()
    {
        if (!_parsed)
        {
            string r = raw;
            if (r !is null && r.length > 0)
            {
                if (source == "pubmed" || source == "pmc")
                    _doc = parseXml(r);
                else
                    _doc = parseWikitext(r);

                _doc.drop(Document.dropFlag(
                    NodeType.Template,
                    NodeType.Reference,
                    NodeType.Comment,
                    NodeType.Category,
                    NodeType.Image,
                ));
            }
            _parsed = true;
        }
        return _doc;
    }

    /// Full plain-text content of all sections.
    string fulltext()
    {
        Document doc = document();
        if (doc.nodes.length == 0)
            return "";
        return doc.sections().map!(s => doc.extractText(s)).join("\n\n");
    }

    /// Plain text of everything before the first section heading.
    string preamble()
    {
        Document doc = document();
        if (doc.nodes.length == 0)
            return "";
        return doc.preambleText();
    }
}

/// Resolve pages for a compound from all enabled sources.
Page[] resolvePage(Compound compound)
{
    Page[] ret;

    if (config.enabledSources.canFind("wikipedia"))
    {
        import akashi.wikipedia : resolvePage_ = resolvePage;
        Page page = resolvePage_(compound);
        if (page !is null)
            ret ~= page;
    }

    if (config.enabledSources.canFind("psychonaut"))
    {
        import akashi.psychonaut : resolvePage_ = resolvePage;
        Page page = resolvePage_(compound);
        if (page !is null)
            ret ~= page;
    }

    return ret;
}
