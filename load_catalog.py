"""Compatibility wrapper for migrated module."""

import warnings

# DEPRECATION SHIM WARNING
warnings.warn(
    "load_catalog.py is a deprecated compatibility shim. Use `workflows.catalog.load_catalog` directly.",
    FutureWarning,
    stacklevel=2,
)


if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.catalog.load_catalog", run_name="__main__")
else:
    from workflows.catalog.load_catalog import *  # noqa: F401,F403
