module tests.integration.psychonaut;

version(psychonaut)
{
    import akashi.page : Page;
    import akashi.psychonaut : Dosage, DosageResult, getDosage, getPages, getPagesByTitle, getReports,
        parseDosageText;
    import akashi.pubchem.compound : Compound, Properties;
    import akashi.text.ast : Document;

    import std.string : indexOf;

private:
    Page[] pages;
    Page[] reports;
    Page page;
    DosageResult dosage;

    static this()
    {
        pages = getPages!"psychonaut"("caffeine", 3);
        assert(pages.length > 0 && pages.length <= 3);
        page = pages[0];
        assert(page.raw.length > 0);
        reports = getReports(page, 5);

        Compound compound = Compound.getOrCreate(2519);
        compound.properties = Properties(
            "Caffeine",
            null,
            null,
            null,
            null,
            194.19,
            double.nan,
            0,
            double.nan,
            -0.1
        );
        dosage = getDosage(compound);
        assert(dosage.dosages.length > 0);
    }

    unittest
    {
        assert(page.source == "psychonaut");
        assert(page.title.length > 0);
    }

    unittest
    {
        assert(page.url.indexOf("https://psychonautwiki.org/wiki/") == 0);
    }

    unittest
    {
        assert(page.raw.length > 500);
    }

    unittest
    {
        Document document = page.document;
        assert(document.nodes.length > 1);
    }

    unittest
    {
        assert(getPagesByTitle!"psychonaut"().length == 0);
    }

    unittest
    {
        assert(reports.length <= 5);
        foreach (report; reports)
            assert(report.source == "psychonaut");
    }

    unittest
    {
        Dosage[] dosages = parseDosageText(
            "| OralROA_Threshold = 10 mg\n"
            ~"| OralROA_Common = [[Dose::50–100 mg]]\n"
        );
        assert(dosages.length == 1);
        assert(dosages[0].route == "Oral");
        assert(dosages[0].threshold == "10 mg");
        assert(dosages[0].common == "50–100 mg");
    }

    unittest
    {
        Dosage[] dosages = parseDosageText(
            "| OralROA_Common = 2 x 10 mg\n"
            ~"| InsufflatedROA_Light = 5 mg<ref>citation</ref>\n"
        );
        assert(dosages.length == 1);
        assert(dosages[0].route == "Insufflated");
        assert(dosages[0].light == "5 mgcitation");
    }

    unittest
    {
        bool foundOral;
        foreach (entry; dosage.dosages)
        {
            if (entry.route == "Oral")
                foundOral = true;
        }
        assert(foundOral);
    }

    unittest
    {
        foreach (entry; dosage.dosages)
            assert(entry.route.length > 0);
    }
}
