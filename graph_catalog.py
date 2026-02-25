"""Compatibility wrapper for migrated module."""

if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.graph_inference.graph_catalog", run_name="__main__")
else:
    from workflows.graph_inference.graph_catalog import *  # noqa: F401,F403
