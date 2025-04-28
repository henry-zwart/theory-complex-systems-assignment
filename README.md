# Theory of Complex Systems assignment

Git repository for Theory of Complex Systems assignment submission.


## Getting started

Clone this repository and its submodules:

```bash
git clone --recurse-submodules git@github.com:henry-zwart/theory-complex-systems-assignment.git
```

Ensure all necessary dependencies are installed:
- [uv](https://docs.astral.sh/uv/): Python package/version management
- [Typst](https://typst.app): Used to write the report

Generate all results, and the corresponing report (solutions file):
```bash
make
```

All figures will be saved to the `results/figures` directory. The compiled report is saved to `solutions.pdf`.

This repository also has a GitHub Action configured to reproduce the results and associated report on any pushes to the main branch. The most recent copy of the solutions pdf can be found under the [GitHub Actions page](https://github.com/henry-zwart/theory-complex-systems-assignment/actions).

## License
This report is licensed under the [MIT License](LICENSE.md).




