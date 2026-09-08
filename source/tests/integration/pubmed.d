module tests.integration.pubmed;

version(pubmed)
{
    import akashi.entrez.eutils : elink;
    import akashi.entrez.pubmed : getPages, getPagesByID;
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
        pages = getPages!"pubmed"("Molegro Virtual Docker", 1);
        assert(pages.length == 1);
        page = pages[0];
        assert(page.raw.length > 0);
        links = elink!("pubmed", "pmc")("31452104");
        assert(links.type == JSONType.object);
    }

    unittest
    {
        assert(pages.length == 1);
        assert(page.source == "pubmed");
    }

    unittest
    {
        assert(page.title.length > 0);
    }

    unittest
    {
        assert(page.url.indexOf("https://pubmed.ncbi.nlm.nih.gov/") == 0);
    }

    unittest
    {
        assert(page.raw.indexOf("<PubmedArticle") >= 0);
        assert(page.raw.indexOf("<Abstract") >= 0);
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
        assert(page.document.source == page.raw);
    }

    unittest
    {
        assert(getPagesByID!"pubmed"().length == 0);
    }

    unittest
    {
        assert("linksets" in links);
    }
}
