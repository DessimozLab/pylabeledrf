# pylabeledrf

pylabeledrf is a heuristic to compute an extension of the Robinson Foulds tree topological distance between labeled topologies, where inner nodes are labeled with speciation or duplication events.

## Citation
If you use our package in your work, please consider citing:

Samuel Briand, Christophe Dessimoz, Nadia El-Mabrouk, Manuel Lafond, Gabriela Lobinska, *Extending the Robinson-Foulds distance to reconciled trees*, submitted


## Installation

The package requires Python 3 (>=3.6). The easiest way to install is using 
`pip`, to install the <a href="https://pypi.org/project/pylabeledrf/">package from 
PyPI</a>.

```
pip install pylabeledrf
```

## Documentation

Documentation is available <a href="http://dessimozlab.github.io/pylabeledrf/build/html/">here</a>.

## Example

```python
from pylabeledrf.computeLRF import *
import dendropy
taxa = dendropy.TaxonNamespace()

# retrieve the test TP53 reconciled tree (from Ensembl compara 96)
p53 = dendropy.Tree.get_from_url(
    'https://raw.githubusercontent.com/DessimozLab/pylabeledrf/master/test/p53.nhx', 
    'newick', taxon_namespace=taxa)
t1 = parseEnsemblLabels(p53)

# introduce 5 random edits and compute the distance
t2 = mutateLabeledTree(t1, 5)
computeLRF(t1,t2)

# randomise the labels and compute the distance
t3 = randomLabels(t1)
computeLRF(t1,t3)
```


## License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
