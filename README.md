# pygeomatch

Python library for matching geospatial vector datasets.
For now, we only have surface matching algorithms but others will come soon.

# installation

## from source

```shell
git clone git@github.com:umrlastig/pymatch.git
cd pymatch
```

[Install uv](https://docs.astral.sh/uv/getting-started/installation/#installation-methods)

Install dependencies:
```shell
uv sync
```

```shell
source .venv/bin/activate
```

# usage

You can now use the *pygeomatch* script to match your data.
You can see the different parameters with:
```shell
pygeomatch --help
```

For instance, with the test data in directory *data* using the *MCA* algorithm and saving the results in *test.gpkg*, you can use the following:
```shell
pygeomatch data/popRef.shp data/popComp.shp test.gpkg MCA
```

# Run tests

```shell
pytest -v --cov  --cov-report term-missing
```

## References

Du, H. (2015). Matching disparate geospatial datasets and validating matches using spatial logic (Doctoral dissertation, University of Nottingham).

Yan, Y., Sun, Y., Wang, S., Lu, Y., Hu, Y., & Lu, M. (2025). Research on Multi-Scale Vector Road-Matching Model Based on ISOD Descriptor. ISPRS International Journal of Geo-Information, 14(7), 280.

## Data

Ground-truth OSM/BDTopo centre Strasbourg (stage) [-> link BiblioDocs]

