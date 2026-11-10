# Contributing

Keep changes focused. Explain what changed and why it was needed.

## Before Making Changes

Read the relevant implementation and tests first. Match the style of neighboring D code.

- Retrieval and request plumbing lives in `source/akashi/request.d`.
- Source-specific APIs live under `source/akashi/pubchem/`, `source/akashi/wikipedia.d`, `source/akashi/psychonaut.d`, and `source/akashi/entrez/`.
- Wikitext and XML parsing lives under `source/akashi/text/`.
- Live integration tests live under `source/tests/integration/`.

## Tests

Run the library build and offline tests before submitting a change.

```sh
dub build --compiler=ldc2 --config=library
dub test --compiler=ldc2 --config=library
```

Add or update tests for behavior changes. Keep live service requests behind their existing version flags; see [Testing](TESTING.md) for the flags and what each suite hits.

## Documentation

Update the README when public behavior changes. Check examples against the code, including argument order and defaults.

## Reporting Issues

Include a minimal example, expected and actual behavior, compiler and DUB versions, and the source involved (PubChem, Wikipedia, PsychonautWiki, PubMed, or PMC). Sanitize logs and never commit credentials, tokens, or response bodies containing private data.
