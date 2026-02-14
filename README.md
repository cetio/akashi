# Akashi
 
Akashi is a D library for fetching and parsing compound-centric content from public knowledge sources.
It provides:
 
- PubChem compound lookups (properties, identifiers, 3D conformers, similarity search)
- Wikipedia and PsychonautWiki page resolution
- PubMed / PMC retrieval (via NCBI Entrez)
- Wikitext + XML parsing into a simple AST for rendering or text extraction
 
Akashi is designed to be used as a standalone dependency in D projects (CLI tools, GUI apps, services).
It is the data and parsing layer used by Chemica.
 
## Install
 
### DUB (recommended)
 
Add Akashi to your `dub.json`:
 
```json
{
	"dependencies": {
		"akashi": "*"
	}
}
```
 
For local development as a path dependency:
 
```json
{
	"dependencies": {
		"akashi": {
			"path": "../akashi"
		}
	}
}
```
 
## Usage
 
### PubChem compound lookup
 
```d
import akashi.pubchem : getProperties;
import akashi.pubchem.compound : Compound;

void main()
{
    Compound c = getProperties("caffeine");
    if (c is null)
        return;

    // Example fields depend on what PubChem returns for the query.
    // Use the Compound you get back directly.
}
```
 
### Resolve pages from enabled sources
 
```d
import akashi.page : resolvePage;
import akashi.pubchem : getProperties;
import akashi.page : Page;
import akashi.pubchem.compound : Compound;

void main()
{
    Compound c = getProperties("caffeine");
    if (c is null)
        return;

    Page[] pages = resolvePage(c);
    foreach (p; pages)
    {
        // p.title, p.source, p.url
        // p.fulltext() for a plain-text view
    }
}
```
 
## Modules
 
- `akashi.pubchem`
  - PubChem compound properties and 3D conformer retrieval.
- `akashi.wikipedia`
  - Wikipedia resolution and content retrieval.
- `akashi.psychonaut`
  - PsychonautWiki access (dosage, reports, pages).
- `akashi.entrez`
  - NCBI Entrez (PubMed + PMC) access.
- `akashi.page`
  - Common `Page` type and page resolution across sources.
- `akashi.text`
  - Wikitext + XML parsing into a small AST (`akashi.text.ast`).
 
## License
 
Akashi is licensed under the [AGPL-3.0 license](LICENSE.txt).
