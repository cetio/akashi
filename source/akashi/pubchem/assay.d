module akashi.pubchem.assay;

import akashi.pubchem : getAssay, getCompoundsByAssay;
import akashi.pubchem.compound;
import akashi.pubchem.protein;

struct AssayDate
{
    int year;
    int month;
    int day;
}

struct AssayCounts
{
    int all;
    int active;
    int inactive;
    int inconclusive;
    int unspecified;
    int probe;
}

struct AssayTarget
{
    string accession;
    string name;

    Protein protein()
    {
        return accession.length == 0 ? null : Protein.getOrCreate(accession);
    }
}

struct AssayResult
{
    Assay assay;
    Compound compound;

    int sid;
    int panelMemberID;
    string outcome;

    string targetAccession;
    int targetGeneID;

    double activityValue = double.nan;
    string activityValueUnit;
    string activityName;

    string assayName;
    string assayType;
    int pubmedID;
    bool rnai;

    string[string] extras;

    bool isActive() const
        => outcome == "Active";

    Protein protein()
    {
        if (targetAccession.length > 0)
            return Protein.getOrCreate(targetAccession);

        if (assay !is null)
        {
            AssayTarget[] assayTargets = assay.targets();
            if (assayTargets.length == 1)
                return assayTargets[0].protein();
        }

        return null;
    }
}

class Assay
{
    static Assay[int] registry;

    static Assay getOrCreate(int aid)
    {
        assert(aid != 0, "AID must not be zero");
        if (auto p = aid in registry)
            return *p;

        Assay assay = new Assay(aid);
        registry[aid] = assay;
        return assay;
    }

package:
    this(int aid)
    {
        assert(aid != 0, "AID must not be zero");
        this.aid = aid;
    }

    string _name;
    string[] _description;
    string[] _protocol;
    string[] _comment;
    AssayTarget[] _targets;
    bool _summaryLoaded;
    Compound[][string] _compounds;
    bool[string] _compoundsLoaded;

public:
    int aid;
    string sourceName;
    string sourceID;
    int numberOfTIDs;
    bool hasScore;
    string method;
    int assayVersion;
    int revision;
    AssayDate lastDataChange;
    AssayCounts sidCounts;
    AssayCounts cidCounts;

    static Assay byID(int aid)
    {
        return getAssay(aid);
    }

    ref string name()
    {
        ensureSummary();
        return _name;
    }

    string[] description()
    {
        ensureSummary();
        return _description;
    }

    string[] protocol()
    {
        ensureSummary();
        return _protocol;
    }

    string[] comment()
    {
        ensureSummary();
        return _comment;
    }

    AssayTarget[] targets()
    {
        ensureSummary();
        return _targets;
    }

    Protein[] proteins()
    {
        Protein[] ret;
        bool[string] seen;

        foreach (target; targets())
        {
            if (target.accession.length == 0 || target.accession in seen)
                continue;

            seen[target.accession] = true;
            ret ~= target.protein();
        }

        return ret;
    }

    Compound[] compounds(string cidsType = "all")
    {
        if (auto p = cidsType in _compoundsLoaded)
        {
            if (*p)
            {
                if (auto compounds = cidsType in _compounds)
                    return *compounds;
                return null;
            }
        }

        _compounds[cidsType] = getCompoundsByAssay(aid, cidsType);
        _compoundsLoaded[cidsType] = true;
        return _compounds[cidsType];
    }

private:
    void ensureSummary()
    {
        if (!_summaryLoaded)
            getAssay(aid);
    }
}
