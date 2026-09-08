module tests.integration.pmc;

version(pmc)
{
    import akashi.entrez.eutils : elink;
    import akashi.entrez.pmc : getDOI, getPages, getPagesByID;
    import akashi.page : Page;
    import akashi.text.ast : Document;

    import std.json : JSONType, JSONValue;
    import std.string : indexOf;

private:
    Page[] pages;
    Page page;
    JSONValue links;

    static this()
    {
        pages = getPages!"pmc"("GenBank", 1);
        assert(pages.length == 1);
        page = pages[0];
        assert(page.raw.length > 0);
        links = elink!("pmc", "pubmed")("3531190");
        assert(links.type == JSONType.object);
    }

    unittest
    {
        assert(pages.length == 1);
        assert(page.source == "pmc");
    }

    unittest
    {
        assert(page.title.length > 0);
    }

    unittest
    {
        assert(page.url.indexOf("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC") == 0);
    }

    unittest
    {
        assert(page.raw.indexOf("<article") >= 0);
        assert(page.raw.indexOf("<article-title") >= 0);
    }

    unittest
    {
        assert(getDOI(page).length > 0);
    }

    unittest
    {
        Document document = page.document;
        assert(document.nodes.length > 1);
    }

    unittest
    {
        Document document = page.document;
        assert(document.sections.length > 0);
    }

    unittest
    {
        assert(page.fulltext.length > 0);
    }

    unittest
    {
        assert(getPagesByID!"pmc"().length == 0);
        assert(getDOI(null) is null);
    }

    unittest
    {
        assert("linksets" in links);
    }
}
