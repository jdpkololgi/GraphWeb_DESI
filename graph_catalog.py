"""Compatibility wrapper for migrated module."""

import warnings

# DEPRECATION SHIM WARNING
warnings.warn(
    "graph_catalog.py is a deprecated compatibility shim. Use `workflows.graph_inference.graph_catalog` directly.",
    FutureWarning,
    stacklevel=2,
)


if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.graph_inference.graph_catalog", run_name="__main__")
else:
    from workflows.graph_inference.graph_catalog import *  # noqa: F401,F403
