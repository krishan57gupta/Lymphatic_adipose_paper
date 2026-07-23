# Path configuration

Several original scripts contain the placeholder `/yourDataAndCodeFolder/` or `yourDataAndCodeFolder`. Before running, replace it with the absolute path to this repository, or refactor the scripts to read `FOXC2_PROJECT_ROOT` from the environment.

Example:

```bash
export FOXC2_PROJECT_ROOT=/absolute/path/to/Lymphatic_adipose_paper
```

The publication release should not contain private user-specific paths.
