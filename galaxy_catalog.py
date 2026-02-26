"""Compatibility wrapper for migrated module."""

import warnings

# DEPRECATION SHIM WARNING
warnings.warn(
    "galaxy_catalog.py is a deprecated compatibility shim. Use `workflows.utilities.galaxy_catalog` directly.",
    FutureWarning,
    stacklevel=2,
)


if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.utilities.galaxy_catalog", run_name="__main__")
else:
    from workflows.utilities.galaxy_catalog import *  # noqa: F401,F403
