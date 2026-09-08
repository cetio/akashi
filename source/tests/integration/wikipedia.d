module tests.integration.wikipedia;

version(wikipedia)
{
    import akashi.page : Page;
    import akashi.pubchem.compound : Compound, Properties;
    import akashi.text.ast : Document;
    import akashi.wikipedia : extractCIDs, getPages, getPagesByTitle, resolvePage;

    import std.algorithm : canFind;
    import std.string : indexOf;

private:
    Page[] pages;
    Page[] titlePages;
    Page page;
    Page resolvedPage;

    static this()
    {
        pages = getPages!"wikipedia"("caffeine", 3);
        assert(pages.length > 0 && pages.length <= 3);
        page = pages[0];
        assert(page.raw.length > 0);
        titlePages = getPagesByTitle!"wikipedia"("Caffeine", "Aspirin");
        assert(titlePages.length == 2);

        Compound compound = Compound.getOrCreate(2519);
        compound.properties = Properties("Caffeine");
        resolvedPage = resolvePage(compound);
        assert(resolvedPage !is null);
    }

    unittest
    {
        assert(page.source == "wikipedia");
        assert(page.title.length > 0);
    }

    unittest
    {
        assert(page.url.indexOf("https://en.wikipedia.org/wiki/") == 0);
    }

    unittest
    {
        assert(page.raw.length > 1_000);
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
        assert(page.preamble.length > 0);
        assert(page.fulltext.length > page.preamble.length);
    }

    unittest
    {
        assert(titlePages.length == 2);
        assert(titlePages[0].source == "wikipedia");
        assert(titlePages[1].source == "wikipedia");
    }

    unittest
    {
        assert(getPagesByTitle!"wikipedia"().length == 0);
    }

    unittest
    {
        Page identifiers = Page.fromRaw(
            "Identifiers",
            "wikipedia",
            "{{Chembox|PubChem = 2519|PubChemCID = 2244}} {{PubChem|3672}} "
                ~"https://pubchem.ncbi.nlm.nih.gov/compound/962"
        );
        string[] cids = extractCIDs(identifiers);
        assert(cids.canFind("2519"));
        assert(cids.canFind("2244"));
        assert(cids.canFind("3672"));
        assert(cids.canFind("962"));
    }

    unittest
    {
        assert(resolvedPage.source == "wikipedia");
        assert(resolvedPage.title.length > 0);
    }
}
