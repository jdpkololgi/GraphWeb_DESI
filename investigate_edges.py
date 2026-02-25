"""Compatibility wrapper for migrated module."""

if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.utilities.investigate_edges", run_name="__main__")
else:
    from workflows.utilities.investigate_edges import *  # noqa: F401,F403
