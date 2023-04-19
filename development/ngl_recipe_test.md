## Recipe

```
rm -rf /opt/conda/share/jupyter/labextensions/nglview-js-widgets
rm -rf /opt/conda/share/jupyter/nbextensions/nglview-js-widgets

pip uninstall ipywidgets
pip uninstall nglview 

conda install ipywidgets==7.6.5 --no-deps
conda install nglview -c conda-forge --no-deps  

jupyter nbextension install widgetsnbextension --py --sys-prefix 
jupyter nbextension enable --py --sys-prefix widgetsnbextension 
jupyter-nbextension enable nglview --py --sys-prefix
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


