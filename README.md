# Repository for actualgroupproject - solving the Time-Dependant Schrodinger Equation

The directory is structured as follows:
```
actualgroupproject/
├── .gitignore
├── LICENSE
├── README.md
├── pyproject.toml
├── requirements.txt
├── setup.cfg
├── .github/
│   └── workflows/
│              └── python_test.yml
├── docs/
│   └── ../
│   └── workflows/
├── src/
│   └── example_package/
│       ├── __init__.py
│       ├── command_line_interface.py
│       └── example.py└── example_package/
│   |__ actualgroupproject_package/
│       ├── command_line_interface.py
│       └── example.py
|       |__solving_tdse.py
└── tests/
        └── test_example.py
|__Group_A_TDSEsol.py

```
Group_A_TDSEsol, when istalled, will allow the user to call the function solve_TDSE() to return an animation of the wavefunction over time for an inputted potential V(x)
