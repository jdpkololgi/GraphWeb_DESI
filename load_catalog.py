"""Compatibility wrapper for migrated module."""

if __name__ == "__main__":
    import runpy

    runpy.run_module("workflows.catalog.load_catalog", run_name="__main__")
else:
    from workflows.catalog.load_catalog import *  # noqa: F401,F403
