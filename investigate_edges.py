"""Compatibility wrapper for migrated module."""

import warnings

# DEPRECATION SHIM WARNING
warnings.warn(
    "investigate_edges.py is a deprecated compatibility shim. Use `workflows.utilities.investigate_edges` directly.",
    FutureWarning,
    stacklevel=2,
)


if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.utilities.investigate_edges", run_name="__main__")
else:
    from workflows.utilities.investigate_edges import *  # noqa: F401,F403
