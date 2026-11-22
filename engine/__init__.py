"""
Graph-theoretic utility network engine (QGIS-independent).

The engine designs water, electricity and gas networks on the plane street
graph of a district and checks every rule of the rule packs as a property of
graphs embedded in the plane (rules/graph_rules.py).  It depends only on
networkx, shapely, geopandas, numpy and scipy, so the QGIS plugin can call it.
"""
