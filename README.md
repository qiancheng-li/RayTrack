# RayTrack

RayTrack contains code, experimental data, and scene assets for radio-map-based
localization and trajectory reconstruction.

## File structure

```text
RayTrack/
|-- assets/
|   |-- raytrack_scene/            Sionna scene XML and building meshes
|   |-- raytrack_model/            Legacy Raytrack scene models
|   `-- uw.osm                     OpenStreetMap source for scene generation
|-- data/
|   |-- algorithm_inputs/          Compact MAT and NPZ algorithm inputs
|   `-- experimental_standardized/ Standardized localization and tracking CSVs
|-- matlab/
|   |-- helpers/                   Coordinate, sampling, peak, and plotting helpers
|   |-- preprocessing/             Localization input-generation scripts
|   |-- tracking/                  Tracking algorithm functions
|   |-- localization_demo.m        MATLAB localization entry point
|   `-- run_tracking_demo.m        MATLAB tracking entry point
|-- python/
|   |-- SceneGenerator.py          Scene and building-mesh generator
|   |-- tracking_algo.py           Python tracking entry point
|   |-- utils.py                   Shared geometry utilities
|   |-- generate_27_frequency_radio_map.ipynb
|   |-- localization.ipynb
|   |-- height_sensitivity.ipynb
|   |-- material_sensitivity.ipynb
|   `-- Sionna_PowerMat_Generator.ipynb
`-- requirements-draft.txt         Python dependencies
```

### Data organization

`data/experimental_standardized/` contains one folder per base station:
`Station-Roof`, `Station-Lab1`, `Station-Lab2`, and `Station-Chem`. Each
station contains ten interpolated tracking trajectories and one combined
fixed-location measurement CSV.

`data/algorithm_inputs/` contains the compact radio-map, coordinate, and recent
path inputs used directly by the retained algorithms.
