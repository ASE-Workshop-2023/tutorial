## Recipe

```
## remove old version installs
rm -rf /opt/conda/share/jupyter/labextensions/nglview-js-widgets
rm -rf /opt/conda/share/jupyter/nbextensions/nglview-js-widgets

## remove pip installs
pip uninstall ipywidgets
pip uninstall nglview 

### ***open up new terminal (or it will not register that packages have been removed)

conda install ipywidgets -c conda-forge --no-deps
conda install -c conda-forge nglview  --no-deps  
```

## Test it

- open notebook

```
from ase import Atoms
import nglview

d = 1.10
molecule = Atoms(['N', 'N'], positions=[(0., 0., 0.), (0., 0., d)])

from ase.visualize import view
view(molecule,viewer='ngl')
```

---> should see a molecule


