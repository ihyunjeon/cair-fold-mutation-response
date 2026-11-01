# src/

Shared code. If two people would otherwise write it twice, it belongs here rather than
in a notebook.

| module | what |
|---|---|
| `mutations.py` | Parse variants, build mutant sequences, and the substitution assertion |
| `geometry.py` | Superpose two structures and measure per-residue displacement |

Import from a notebook in `notebooks/`:

```python
import sys; sys.path.append('../src')
from mutations import apply_variant, assert_variant_set
```
