#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import gzip
from typing import Iterable, Tuple, Dict, Set, List

def smart_open(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode, encoding=None if "b" in mode else "utf-8")

def _rc(orient: str) -> str:
    return '+' if orient == '-' else '-'

def _split_sid_and_orient(sid: str, orient_field: str = None):
    if orient_field in ('+', '-'):
        return sid, orient_field
    if sid and sid[-1] in ('+', '-'):
        return sid[:-1], sid[-1]
    return sid, None

def parse_gfa_oriented(gfa_path: str, add_rc_edges: bool):

    segments: Set[str] = set()
    nodes: Set[str] = set()
    edges: List[Tuple[str, str]] = []

    with smart_open(gfa_path, "rt") as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n\r").split("\t")
            if not parts:
                continue
            t = parts[0]
            if t == "S":
                # GFA1: S <sid> <seq> ...
                # GFA2: S <sid> <len> <seq> ...
                if len(parts) >= 2:
                    sid = parts[1]
                    segments.add(sid)
                    nodes.add(sid + '+')
                    nodes.add(sid + '-')
            elif t == "L":
                # L <from> <from_orient> <to> <to_orient> <overlap> ...
                if len(parts) >= 6:
                    u_raw, o1_raw = parts[1], parts[2]
                    v_raw, o2_raw = parts[3], parts[4]
                    u_sid, o1 = _split_sid_and_orient(u_raw, o1_raw)
                    v_sid, o2 = _split_sid_and_orient(v_raw, o2_raw)
                    if o1 in ('+', '-') and o2 in ('+', '-'):
                        src = u_sid + o1
                        dst = v_sid + o2
                        edges.append((src, dst))
                        if add_rc_edges:
                            edges.append((v_sid + _rc(o2), u_sid + _rc(o1)))
                        nodes.add(u_sid + '+'); nodes.add(u_sid + '-')
                        nodes.add(v_sid + '+'); nodes.add(v_sid + '-')
            elif t == "E":
                # GFA2 heuristic: looking for s1,o1,s2,o2
                s1 = o1 = s2 = o2 = None
                if len(parts) >= 6 and parts[3] in ('+', '-') and parts[5] in ('+', '-'):
                    s1_raw, o1_raw = parts[2], parts[3]
                    s2_raw, o2_raw = parts[4], parts[5]
                    s1, o1 = _split_sid_and_orient(s1_raw, o1_raw)
                    s2, o2 = _split_sid_and_orient(s2_raw, o2_raw)
                elif len(parts) >= 4:
                    s1, o1 = _split_sid_and_orient(parts[2], None)
                    s2, o2 = _split_sid_and_orient(parts[3], None)
                if s1 and s2 and o1 in ('+', '-') and o2 in ('+', '-'):
                    src = s1 + o1
                    dst = s2 + o2
                    edges.append((src, dst))
                    if add_rc_edges:
                        edges.append((s2 + _rc(o2), s1 + _rc(o1)))
                    nodes.add(s1 + '+'); nodes.add(s1 + '-')
                    nodes.add(s2 + '+'); nodes.add(s2 + '-')
    return segments, nodes, edges

def build_index(nodes: Set[str], edges: Iterable[Tuple[str, str]], include_isolated: bool) -> Dict[str, int]:
    """
    Builds a mapping from oriented_node_id -> integer index (0..n-1).
    - include_isolated=True: includes all known oriented nodes
    - otherwise: includes only those present in at least one edge
    """

    if include_isolated:
        node_set: Set[str] = set(nodes)
    else:
        node_set = set()
        for u, v in edges:
            node_set.add(u); node_set.add(v)
    ordered = sorted(node_set)
    return {nid: i for i, nid in enumerate(ordered)}

def write_sgraph(out_path: str, index: Dict[str, int], edges: Iterable[Tuple[str, str]], one_based: bool, dedup_edges: bool):
    """
    Writes the .sgraph:
    - first line: n m
    - then m lines: a b
    If dedup_edges=True, removes exact multi-edges.
    """

    if dedup_edges:
        conv_edges_set: Set[Tuple[int, int]] = set()
        for u, v in edges:
            if u in index and v in index:
                a = index[u]; b = index[v]
                if one_based:
                    a += 1; b += 1
                conv_edges_set.add((a, b))
        conv_edges = sorted(conv_edges_set)
    else:
        conv_edges: List[Tuple[int, int]] = []
        for u, v in edges:
            if u in index and v in index:
                a = index[u]; b = index[v]
                if one_based:
                    a += 1; b += 1
                conv_edges.append((a, b))
        conv_edges.sort()

    n = len(index)
    m = len(conv_edges)
    with smart_open(out_path, "wt") as out:
        out.write(f"{n} {m}\n")
        for a, b in conv_edges:
            out.write(f"{a} {b}\n")

def write_mapping(map_path: str, index: Dict[str, int], one_based: bool):
    """
    Writes a mapping file: oriented_node<tab>index
    """
    inverse = sorted(index.items(), key=lambda kv: kv[1])
    with smart_open(map_path, "wt") as out:
        for nid, i in inverse:
            if one_based:
                i += 1
            out.write(f"{nid}\t{i}\n")

def main():
    ap = argparse.ArgumentParser(description="Converts a .gfa into an oriented .sgraph (n m followed by m lines a b).")
    ap.add_argument("gfa", help="Path to the .gfa file (possibly .gz)")
    ap.add_argument("sgraph", help="Output path for the .sgraph file (possibly .gz)")
    ap.add_argument("--edges-only", action="store_true",
                    help="Include only oriented nodes present in edges (ignore isolated segments).")
    ap.add_argument("--one-based", action="store_true",
                    help="Use 1..n indexing instead of 0..n-1.")
    ap.add_argument("--map-out", help="Write a mapping from oriented_node (sid+/-) -> index (for reference).")
    ap.add_argument("--add-rc-edges", action="store_true",
                    help="Also add implicit reverse-complement edges (may double bubbles).")
    ap.add_argument("--dedup-edges", action="store_true",
                    help="Remove exact multi-edges in the output.")
    args = ap.parse_args()


    segments, nodes, edges = parse_gfa_oriented(args.gfa, add_rc_edges=args.add_rc_edges)
    dedup_edges = args.dedup_edges
    index = build_index(nodes, edges, include_isolated=not args.edges_only)
    write_sgraph(args.sgraph, index, edges, one_based=args.one_based, dedup_edges=dedup_edges)
    if args.map_out:
        write_mapping(args.map_out, index, one_based=args.one_based)

if __name__ == "__main__":
    main()