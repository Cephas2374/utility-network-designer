"""
Storage of designs and analysis results (model.pkl, water_p2.pkl, elec_p3.pkl, corridor_p4.pkl).

A design holds networkx graphs whose edges carry shapely geometries, GeoDataFrames and the engine's own
StreetGraph and CorridorGraph, which have no faithful plain-text form, so it is stored with pickle.
Unpickling an arbitrary file can run code, so every file is read through RestrictedUnpickler: it
resolves only the classes a design is made of (numpy, pandas, pyarrow, geopandas, shapely, networkx,
pyproj, the engine's own modules and two plain builtins) and refuses everything else (os, subprocess,
builtins.eval, ...).  A tampered design folder therefore raises UnsafeFileError instead of running code.

Added 2026-11-23 for the QGIS plugin repository's security scan: the only place where the engine
imports pickle.
"""
import pickle  # nosec B403 - only used through RestrictedUnpickler below (allow-list of classes)

# modules (and their sub-modules) whose classes a design may contain
ALLOWED_MODULES = ("numpy", "pandas", "pyarrow", "geopandas", "shapely", "networkx", "pyproj", "engine",
                   "collections", "datetime")
# plain builtins that appear in designs (no callables that execute code)
ALLOWED_BUILTINS = {"bytearray", "bytes", "complex", "dict", "frozenset", "list", "range", "set", "slice",
                    "tuple", "float", "int", "str", "bool"}


class UnsafeFileError(pickle.UnpicklingError):
    """The file asks for a class outside the allow-list."""


class RestrictedUnpickler(pickle.Unpickler):  # nosec B301 - find_class enforces the allow-list
    def find_class(self, module, name):
        root = module.split(".")[0]
        if module == "builtins" and name in ALLOWED_BUILTINS:
            return super().find_class(module, name)
        if root in ALLOWED_MODULES:
            return super().find_class(module, name)
        raise UnsafeFileError(f"refused to load {module}.{name}: not part of a design")


def load(path):
    """Read a stored design or result through the allow-list."""
    with open(path, "rb") as fh:
        return RestrictedUnpickler(fh).load()


def save(obj, path):
    """Write a design or result."""
    with open(path, "wb") as fh:
        pickle.dump(obj, fh)  # same protocol as before the change
