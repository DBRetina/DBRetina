"""DBRetina viz — interactive Plotly visualizations for pairwise data."""

import os
from collections import Counter
from typing import Optional

import click
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


# ── Helper Functions ──────────────────────────────────────────────


def _build_igraph(graph, metric: str):
    """Return a PairwiseGraph's already-built igraph plus its id maps.

    The graph is built once in ``PairwiseGraph.__init__`` with vertex ``i``
    keyed by ``sorted(names_map)[i]`` and edge ``weight``/``shared_features``
    attributes; callers here treat it read-only.

    Returns:
        (ig_graph, all_ids, id_idx, names_map)
        where all_ids is the sorted list of group IDs,
        id_idx maps group_id -> igraph vertex index,
        names_map maps group_id -> group_name.
    """
    names_map = graph._names_map
    all_ids = graph._idx_to_gid
    id_idx = graph._gid_to_idx
    return graph.to_igraph(), all_ids, id_idx, names_map


def _compute_layout(ig_graph, large_threshold=500):
    """Compute 2D layout coordinates.

    Uses Fruchterman-Reingold for small graphs, DrL for large ones.
    """
    n = ig_graph.vcount()
    has_weights = ig_graph.ecount() > 0 and "weight" in ig_graph.es.attributes()

    if n > large_threshold:
        layout = ig_graph.layout_drl(weights="weight" if has_weights else None)
    else:
        layout = ig_graph.layout_fruchterman_reingold(
            weights="weight" if has_weights else None,
            niter=500,
        )
    return layout


def _community_colors(communities: dict[str, int]) -> tuple[dict[str, str], dict[int, str]]:
    """Map community IDs to distinct Plotly colors.

    Returns:
        (name_to_color, cid_to_color)
    """
    palette = (
        px.colors.qualitative.D3
        + px.colors.qualitative.G10
        + px.colors.qualitative.T10
        + px.colors.qualitative.Alphabet
    )
    unique_ids = sorted(set(communities.values()))
    cid_to_color = {cid: palette[i % len(palette)] for i, cid in enumerate(unique_ids)}
    name_to_color = {name: cid_to_color[cid] for name, cid in communities.items()}
    return name_to_color, cid_to_color


def _scale_sizes(values: dict, min_size: float = 5, max_size: float = 25) -> dict:
    """Linear-scale a dict of float values to a pixel size range."""
    if not values:
        return {}
    vals = list(values.values())
    lo, hi = min(vals), max(vals)
    if hi == lo:
        return {k: (min_size + max_size) / 2 for k in values}
    return {
        k: min_size + (v - lo) / (hi - lo) * (max_size - min_size)
        for k, v in values.items()
    }


def _generate_index_html(files: list[tuple[str, str]], output_dir: str):
    """Write an index.html that links to all generated figures.

    Args:
        files: List of (filename, description) tuples.
        output_dir: Directory to write index.html into.
    """
    rows = "\n".join(
        f'<li><a href="{fname}">{desc}</a></li>' for fname, desc in files
    )
    html = f"""<!DOCTYPE html>
<html><head><title>DBRetina Visualizations</title>
<style>
body {{ font-family: sans-serif; max-width: 600px; margin: 40px auto; }}
a {{ color: #1f77b4; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
li {{ margin: 8px 0; }}
</style></head>
<body>
<h1>DBRetina Visualizations</h1>
<ul>{rows}</ul>
</body></html>"""
    with open(os.path.join(output_dir, "index.html"), "w") as f:
        f.write(html)


# ── Figure 1: Pairwise Summary ───────────────────────────────────


def plot_pairwise_summary(store, metric: str = "ochiai") -> go.Figure:
    """Dataset overview: metric distributions, top pairs, degree distribution.

    Args:
        store: PairwiseStore instance.
        metric: Metric to highlight for top pairs.

    Returns:
        Plotly Figure with 2x2 subplots.
    """
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            "Metric Distributions",
            f"Top 20 Pairs by {metric}",
            "Degree Distribution",
            "Dataset Overview",
        ],
        vertical_spacing=0.12,
        horizontal_spacing=0.1,
    )

    # Subplot 1: Violin per metric (sampled)
    metrics = ("containment", "ochiai", "jaccard", "csi", "dice")
    for m in metrics:
        try:
            sample_df = store.sql(
                f"SELECT {m} FROM pairs USING SAMPLE 5000"
            ).fetchdf()
            fig.add_trace(
                go.Violin(
                    y=sample_df[m], name=m, box_visible=True,
                    meanline_visible=True, showlegend=False,
                ),
                row=1, col=1,
            )
        except Exception:
            pass

    # Subplot 2: Top-20 horizontal bar
    top_df = store.top_pairs(metric, n=20)
    if len(top_df) > 0:
        labels = []
        for _, r in top_df.iterrows():
            n1 = str(r.get("group_1_name", r.get("group_1_id", "")))
            n2 = str(r.get("group_2_name", r.get("group_2_id", "")))
            labels.append(f"{n1[:25]} / {n2[:25]}")
        fig.add_trace(
            go.Bar(
                x=top_df[metric].tolist(),
                y=labels,
                orientation="h",
                marker_color="#1f77b4",
                showlegend=False,
            ),
            row=1, col=2,
        )
    fig.update_yaxes(autorange="reversed", row=1, col=2)

    # Subplot 3: Degree distribution (log-log)
    degree_counts = store.group_pair_counts()
    if degree_counts:
        deg_hist = Counter(degree_counts.values())
        degs = sorted(deg_hist.keys())
        counts = [deg_hist[d] for d in degs]
        fig.add_trace(
            go.Scatter(
                x=degs, y=counts, mode="markers",
                marker=dict(size=4, color="#ff7f0e"),
                showlegend=False,
            ),
            row=2, col=1,
        )
        fig.update_xaxes(type="log", title_text="Degree", row=2, col=1)
        fig.update_yaxes(type="log", title_text="Count", row=2, col=1)

    # Subplot 4: Stats text
    manifest = store.manifest
    stats_text = (
        f"<b>Groups:</b> {store.num_groups:,}<br>"
        f"<b>Pairs:</b> {store.num_pairs:,}<br>"
        f"<b>Population:</b> {manifest.get('population_size', '?'):,}<br>"
    )
    if manifest.get("cutoff_metric"):
        stats_text += (
            f"<b>Cutoff:</b> {manifest['cutoff_metric']} "
            f">= {manifest.get('cutoff_threshold', 0)}<br>"
        )
    # Add metric summaries
    for m in ("ochiai", "jaccard", "containment"):
        try:
            s = store.metric_summary(m)
            stats_text += (
                f"<b>{m}:</b> mean={s['mean']:.1f}, "
                f"median={s['median']:.1f}, max={s['max']:.1f}<br>"
            )
        except Exception:
            pass

    fig.add_annotation(
        text=stats_text, xref="x4", yref="y4",
        x=0.5, y=0.5, showarrow=False,
        font=dict(size=13), align="left",
        row=2, col=2,
    )
    fig.update_xaxes(visible=False, row=2, col=2)
    fig.update_yaxes(visible=False, row=2, col=2)

    fig.update_layout(
        height=700, width=1100,
        title_text="Pairwise Dataset Summary",
        template="plotly_white",
    )
    return fig


# ── Figure 2: Community Network ──────────────────────────────────


def plot_community_network(
    store,
    metric: str = "ochiai",
    cutoff: float = 20,
    max_nodes: int = 500,
    max_communities: int = 15,
) -> go.Figure:
    """Disease community network with force-directed layout.

    Nodes colored by Leiden community, sized by degree.
    Hover shows disease name, community, degree, PageRank.

    Args:
        store: PairwiseStore instance.
        metric: Similarity metric for edges.
        cutoff: Minimum metric value.
        max_nodes: Maximum nodes before pruning.
        max_communities: Keep only top N communities if too large.

    Returns:
        Plotly Figure.
    """
    from .pairwise_graph import PairwiseGraph

    graph = PairwiseGraph(store, metric, cutoff)
    actual_cutoff = cutoff

    try:
        # First, get only connected nodes (those with at least one edge)
        edges_df = graph.edges_dataframe()
        names_map = graph._names_map
        connected_names = {
            names_map[gid]
            for gid in set(edges_df["src"].tolist() + edges_df["dst"].tolist())
            if gid in names_map
        }

        if len(connected_names) > max_nodes:
            # Run community detection to prune intelligently
            communities = graph.community_detection(group_names=list(connected_names))
            # Filter to only connected nodes
            communities = {n: c for n, c in communities.items() if n in connected_names}
            comm_sizes = Counter(communities.values())

            # Keep top communities until we reach max_nodes
            keep_names = []
            for cid, _ in comm_sizes.most_common():
                members = [n for n, c in communities.items() if c == cid]
                if len(keep_names) + len(members) > max_nodes and keep_names:
                    break
                keep_names.extend(members)
                if len(keep_names) >= max_nodes:
                    break

            if not keep_names:
                keep_names = list(connected_names)[:max_nodes]

            sub = graph.subgraph(keep_names)
            graph.close()
            graph = sub
        elif len(connected_names) < graph.num_nodes:
            # Remove isolated nodes by creating subgraph of connected ones
            sub = graph.subgraph(list(connected_names))
            graph.close()
            graph = sub

        communities = graph.community_detection()
        pr = graph.pagerank()

        ig_graph, all_ids, id_idx, names_map = _build_igraph(graph, metric)
        layout = _compute_layout(ig_graph, large_threshold=max_nodes)

        degrees = ig_graph.degree()
        names = [names_map[gid] for gid in all_ids]

        name_colors, cid_colors = _community_colors(communities)
        degree_dict = {n: degrees[id_idx[gid]] for gid, n in zip(all_ids, names)}
        node_sizes = _scale_sizes(degree_dict, 5, 25)

        # Edge trace
        edge_x, edge_y = [], []
        for e in ig_graph.es:
            x0, y0 = layout[e.source]
            x1, y1 = layout[e.target]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            mode="lines",
            line=dict(width=0.3, color="rgba(150,150,150,0.3)"),
            hoverinfo="none",
            showlegend=False,
        )

        # Node trace
        node_x = [layout[id_idx[gid]][0] for gid in all_ids]
        node_y = [layout[id_idx[gid]][1] for gid in all_ids]
        node_color = [name_colors.get(n, "#999") for n in names]
        node_size = [node_sizes.get(n, 8) for n in names]
        hover_texts = [
            f"<b>{n}</b><br>"
            f"Community: {communities.get(n, '?')}<br>"
            f"Degree: {degree_dict.get(n, 0)}<br>"
            f"PageRank: {pr.get(n, 0):.4f}"
            for n in names
        ]

        use_gl = len(names) > 300
        scatter_cls = go.Scattergl if use_gl else go.Scatter

        node_trace = scatter_cls(
            x=node_x, y=node_y,
            mode="markers+text",
            marker=dict(size=node_size, color=node_color, line=dict(width=0.5, color="white")),
            text=[n[:20] if degree_dict.get(n, 0) > 5 else "" for n in names],
            textposition="top center",
            textfont=dict(size=7),
            hovertext=hover_texts,
            hoverinfo="text",
            showlegend=False,
        )

        # Community legend traces
        legend_traces = []
        comm_sizes = Counter(communities.values())
        for cid, count in comm_sizes.most_common(max_communities):
            legend_traces.append(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=8, color=cid_colors[cid]),
                name=f"Community {cid} ({count})",
                showlegend=True,
            ))

        fig = go.Figure(data=[edge_trace, node_trace] + legend_traces)
        fig.update_layout(
            title=f"Disease Community Network ({metric} >= {actual_cutoff}, "
                  f"{len(names)} nodes, {ig_graph.ecount()} edges)",
            showlegend=True,
            legend=dict(x=1.02, y=1, font=dict(size=10)),
            hovermode="closest",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            template="plotly_white",
            height=800, width=1100,
            margin=dict(l=20, r=200, t=60, b=20),
        )
        return fig
    finally:
        graph.close()


# ── Figure 3: Neighborhood Spotlight ─────────────────────────────


def plot_neighborhood(
    store,
    group_name: str,
    metric: str = "ochiai",
    cutoff: float = 20,
    hops: int = 2,
) -> go.Figure:
    """Neighborhood spotlight with gene overlay.

    Target disease centered, N-hop neighbors shown.
    Node size = PageRank, color = community, opacity fades by hop distance.
    Edge hover shows shared genes (if dbri_path is set on store).

    Args:
        store: PairwiseStore instance.
        group_name: Target disease name.
        metric: Similarity metric for edges.
        cutoff: Minimum metric value.
        hops: Number of hops from target.

    Returns:
        Plotly Figure.
    """
    from .pairwise_graph import PairwiseGraph

    graph = PairwiseGraph(store, metric, cutoff)
    sub = graph.subgraph_around(group_name, hops=hops)

    try:
        pr = sub.pagerank()
        communities = sub.community_detection()

        ig_graph, all_ids, id_idx, names_map = _build_igraph(sub, metric)
        layout = _compute_layout(ig_graph)

        names = [names_map[gid] for gid in all_ids]

        # Compute hop distances from center via igraph BFS
        target_key = group_name.lower()
        center_gid = None
        for gid in all_ids:
            if names_map[gid].lower() == target_key:
                center_gid = gid
                break

        hop_distances = {}
        if center_gid is not None and center_gid in id_idx:
            center_idx = id_idx[center_gid]
            dists = ig_graph.shortest_paths(source=center_idx)[0]
            for i, gid in enumerate(all_ids):
                d = dists[i]
                hop_distances[names_map[gid]] = d if d != float("inf") else hops + 1
        else:
            hop_distances = {n: 0 for n in names}

        # Gene overlay: shared genes per edge
        has_genes = store._dbri_path is not None
        edge_gene_text = {}
        if has_genes:
            try:
                store._load_gene_sets()
                for e in ig_graph.es:
                    src_name = names_map[all_ids[e.source]]
                    dst_name = names_map[all_ids[e.target]]
                    try:
                        shared = store.shared_features(src_name, dst_name)
                        top5 = sorted(shared)[:5]
                        edge_gene_text[(e.source, e.target)] = ", ".join(top5)
                    except Exception:
                        edge_gene_text[(e.source, e.target)] = ""
            except Exception:
                has_genes = False

        # Gene overlay: top genes per node from raw gene sets
        node_gene_text = {}
        if has_genes:
            try:
                gene_sets = store._load_gene_sets()
                for n in names:
                    genes = gene_sets.get(n.lower(), [])[:5]
                    node_gene_text[n] = ", ".join(genes) if genes else ""
            except Exception:
                pass

        # Colors and sizes
        name_colors, cid_colors = _community_colors(communities)
        pr_sizes = _scale_sizes(pr, 8, 30)

        # Opacity by hop distance
        def _hop_opacity(d):
            if d == 0:
                return 1.0
            elif d == 1:
                return 0.85
            elif d == 2:
                return 0.55
            return 0.35

        # Edge traces: line trace + midpoint hover trace
        edge_x, edge_y = [], []
        mid_x, mid_y, mid_hover = [], [], []

        for e in ig_graph.es:
            x0, y0 = layout[e.source]
            x1, y1 = layout[e.target]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

            # Midpoint for hover
            mx, my = (x0 + x1) / 2, (y0 + y1) / 2
            mid_x.append(mx)
            mid_y.append(my)

            src_name = names_map[all_ids[e.source]]
            dst_name = names_map[all_ids[e.target]]
            weight = e["weight"] if "weight" in e.attributes() else 0
            genes = edge_gene_text.get((e.source, e.target), "")

            hover = f"<b>{src_name}</b> ↔ <b>{dst_name}</b><br>{metric}: {weight:.1f}"
            if genes:
                hover += f"<br>Shared genes: {genes}"
            mid_hover.append(hover)

        edge_line_trace = go.Scatter(
            x=edge_x, y=edge_y,
            mode="lines",
            line=dict(width=0.8, color="rgba(150,150,150,0.4)"),
            hoverinfo="none",
            showlegend=False,
        )

        edge_hover_trace = go.Scatter(
            x=mid_x, y=mid_y,
            mode="markers",
            marker=dict(size=8, color="rgba(0,0,0,0)"),
            hovertext=mid_hover,
            hoverinfo="text",
            showlegend=False,
        )

        # Node trace
        node_x = [layout[id_idx[gid]][0] for gid in all_ids]
        node_y = [layout[id_idx[gid]][1] for gid in all_ids]
        node_color = [name_colors.get(n, "#999") for n in names]
        node_size = [pr_sizes.get(n, 10) for n in names]
        node_opacity = [_hop_opacity(hop_distances.get(n, hops + 1)) for n in names]

        hover_texts = []
        for n in names:
            h = (
                f"<b>{n}</b><br>"
                f"Community: {communities.get(n, '?')}<br>"
                f"Hop: {hop_distances.get(n, '?')}<br>"
                f"PageRank: {pr.get(n, 0):.4f}"
            )
            genes = node_gene_text.get(n, "")
            if genes:
                h += f"<br>Top genes: {genes}"
            hover_texts.append(h)

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode="markers+text",
            marker=dict(
                size=node_size,
                color=node_color,
                opacity=node_opacity,
                line=dict(width=1, color="white"),
            ),
            text=[n[:25] for n in names],
            textposition="top center",
            textfont=dict(size=8),
            hovertext=hover_texts,
            hoverinfo="text",
            showlegend=False,
        )

        # Legend
        legend_traces = []
        comm_sizes = Counter(communities.values())
        for cid, count in comm_sizes.most_common(15):
            legend_traces.append(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=8, color=cid_colors[cid]),
                name=f"Community {cid} ({count})",
            ))

        fig = go.Figure(
            data=[edge_line_trace, edge_hover_trace, node_trace] + legend_traces
        )
        fig.update_layout(
            title=(
                f'Neighborhood of "{group_name}" '
                f"({metric} >= {cutoff}, {hops} hops, "
                f"{len(names)} nodes, {ig_graph.ecount()} edges)"
            ),
            showlegend=True,
            legend=dict(x=1.02, y=1, font=dict(size=10)),
            hovermode="closest",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            template="plotly_white",
            height=800, width=1100,
            margin=dict(l=20, r=200, t=60, b=20),
        )
        return fig
    finally:
        sub.close()
        graph.close()


# ── Figure 4: Threshold Sweep ────────────────────────────────────


def plot_threshold_sweep(
    store,
    group_name: str,
    metric: str = "ochiai",
    cutoff: float = 20,
    thresholds: Optional[list[float]] = None,
    hops: int = 2,
) -> go.Figure:
    """Animated network showing how neighborhood changes with cutoff.

    Uses Plotly slider for interactive threshold control.
    Nodes keep stable positions (anchored from lowest cutoff layout).

    Args:
        store: PairwiseStore instance.
        group_name: Target disease name.
        metric: Similarity metric.
        thresholds: List of cutoff values. Default: [10,20,...,80].
        hops: Neighborhood hops.

    Returns:
        Plotly Figure with animation frames.
    """
    from .pairwise_graph import PairwiseGraph

    if thresholds is None:
        # Start from cutoff (or cutoff-10 to show some growth), step to 80
        start = max(10, int(cutoff) - 10)
        thresholds = list(range(start, 81, 10))
        if not thresholds:
            thresholds = [int(cutoff)]

    thresholds = sorted(thresholds)

    # Build ONE graph at lowest cutoff — extract neighborhood once
    base_graph = PairwiseGraph(store, metric, thresholds[0])
    try:
        sub = base_graph.subgraph_around(group_name, hops=hops)
    except Exception:
        base_graph.close()
        raise

    try:
        # Get all edges with weights from the neighborhood subgraph,
        # mapping group ids to names (downstream filters on the `src`/`dst`
        # name columns and the `weight` column).
        sub_names_map = sub._names_map
        raw_edges = sub.edges_dataframe()
        edges_df = raw_edges.assign(
            src=raw_edges["src"].map(sub_names_map),
            dst=raw_edges["dst"].map(sub_names_map),
        )[["src", "dst", "weight"]]
        all_node_names = list(sub._names_map.values())

        # Compute layout once on the full (lowest cutoff) neighborhood
        ig_full, all_ids, id_idx, names_map = _build_igraph(sub, metric)
        layout = _compute_layout(ig_full)

        # Position map: name -> (x, y)
        pos_map = {}
        for i, gid in enumerate(all_ids):
            pos_map[names_map[gid]] = (layout[i][0], layout[i][1])

        # Community detection on the full neighborhood
        base_communities = sub.community_detection()
        name_colors, _ = _community_colors(base_communities)

    finally:
        sub.close()
        base_graph.close()

    # Build frames by filtering edges in Python
    frames = []
    initial_data = None

    for t in thresholds:
        # Filter edges above this threshold
        if not edges_df.empty:
            filtered = edges_df[edges_df["weight"] >= t]
        else:
            filtered = edges_df

        # Connected nodes at this threshold
        if not filtered.empty:
            connected = set(filtered["src"].tolist() + filtered["dst"].tolist())
        else:
            connected = set()

        # Always include target
        connected.add(group_name)
        frame_names = [n for n in all_node_names if n in connected]

        # Edges using anchor positions
        edge_x, edge_y = [], []
        for _, row in filtered.iterrows():
            src, dst = row["src"], row["dst"]
            if src in pos_map and dst in pos_map:
                x0, y0 = pos_map[src]
                x1, y1 = pos_map[dst]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

        # Node positions
        node_x = [pos_map.get(n, (0, 0))[0] for n in frame_names]
        node_y = [pos_map.get(n, (0, 0))[1] for n in frame_names]
        node_color = [name_colors.get(n, "#999") for n in frame_names]

        n_nodes = len(frame_names)
        n_edges = len(filtered)

        # Simple community count from base communities
        frame_comms = set()
        for n in frame_names:
            if n in base_communities:
                frame_comms.add(base_communities[n])
        n_communities = len(frame_comms)

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y, mode="lines",
            line=dict(width=0.8, color="rgba(150,150,150,0.4)"),
            hoverinfo="none", showlegend=False,
        )
        node_trace = go.Scatter(
            x=node_x, y=node_y, mode="markers+text",
            marker=dict(size=10, color=node_color, line=dict(width=0.5, color="white")),
            text=[n[:20] for n in frame_names],
            textposition="top center", textfont=dict(size=7),
            hovertext=[
                f"<b>{n}</b><br>Community: {base_communities.get(n, '?')}"
                for n in frame_names
            ],
            hoverinfo="text", showlegend=False,
        )

        frame = go.Frame(
            data=[edge_trace, node_trace],
            name=str(int(t)),
            layout=go.Layout(annotations=[dict(
                text=f"{metric} >= {t} | Nodes: {n_nodes} | "
                     f"Edges: {n_edges} | Communities: {n_communities}",
                xref="paper", yref="paper", x=0.5, y=1.05,
                showarrow=False, font=dict(size=14),
            )]),
        )
        frames.append(frame)

        if initial_data is None:
            initial_data = [edge_trace, node_trace]

    if initial_data is None:
        initial_data = [
            go.Scatter(x=[], y=[], mode="lines"),
            go.Scatter(x=[], y=[], mode="markers"),
        ]

    # Slider
    sliders = [dict(
        active=0,
        currentvalue=dict(prefix=f"{metric} >= ", font=dict(size=14)),
        steps=[
            dict(
                label=str(int(t)),
                method="animate",
                args=[[str(int(t))], dict(
                    mode="immediate",
                    frame=dict(duration=300, redraw=True),
                    transition=dict(duration=200),
                )],
            )
            for t in thresholds
        ],
        x=0.1, len=0.8,
        xanchor="left", y=-0.05,
    )]

    fig = go.Figure(data=initial_data, frames=frames)
    fig.update_layout(
        title=f'Threshold Sweep: "{group_name}" ({hops}-hop neighborhood)',
        sliders=sliders,
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        template="plotly_white",
        height=800, width=1000,
        margin=dict(l=20, r=20, t=60, b=80),
    )
    return fig


# ── CLI Command ──────────────────────────────────────────────────


@cli.command(name="viz", epilog=dbretina_doc.doc_url("viz"), help_priority=18)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.option(
    "-i", "--index-prefix", "index_prefix", default=None, type=str,
    help="Index file prefix (for gene overlay on neighborhood)",
)
@click.option(
    "-m", "--metric", default="ochiai", show_default=True, type=str,
    help="Similarity metric for edges",
)
@click.option(
    "-c", "--cutoff", default=20.0, show_default=True, type=float,
    help="Minimum metric value",
)
@click.option(
    "--target", default=None, type=str,
    help="Target disease for neighborhood and threshold sweep",
)
@click.option(
    "--hops", default=2, show_default=True, type=int,
    help="Neighborhood hops",
)
@click.option(
    "--thresholds", "thresholds_str", default=None, type=str,
    help="Comma-separated thresholds for sweep (default: 10,20,...,80)",
)
@click.option(
    "-o", "--output", required=True, type=click.Path(),
    help="Output directory for HTML files",
)
@click.pass_context
def main(ctx, data_path, index_prefix, metric, cutoff, target, hops,
         thresholds_str, output):
    """Generate interactive Plotly visualizations of pairwise data.

    Always generates community_network.html and pairwise_summary.html.
    When --target is given, also generates neighborhood and threshold sweep.

    \b
    Example:
        DBRetina viz -d experiment_pairwise -m ochiai -c 30 -o viz_output/
        DBRetina viz -d experiment_pairwise -i myindex \\
            --target "alzheimer's disease" -o viz_output/
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore

    data_path = os.path.abspath(data_path)
    dbri_path = os.path.abspath(f"{index_prefix}.dbri") if index_prefix else None

    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path, dbri_path=dbri_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data at {data_path}: {e}")
            return
    elif dbri_path:
        store._dbri_path = dbri_path

    os.makedirs(output, exist_ok=True)
    generated_files = []

    try:
        # Always: pairwise summary
        LOGGER.INFO("Generating pairwise summary...")
        fig = plot_pairwise_summary(store, metric=metric)
        path = os.path.join(output, "pairwise_summary.html")
        fig.write_html(path)
        generated_files.append(("pairwise_summary.html", "Pairwise Dataset Summary"))
        LOGGER.INFO(f"  -> {path}")

        # Always: community network
        LOGGER.INFO("Generating community network...")
        fig = plot_community_network(store, metric=metric, cutoff=cutoff)
        path = os.path.join(output, "community_network.html")
        fig.write_html(path)
        generated_files.append(("community_network.html", "Disease Community Network"))
        LOGGER.INFO(f"  -> {path}")

        # If target: neighborhood
        if target:
            LOGGER.INFO(f'Generating neighborhood for "{target}"...')
            fig = plot_neighborhood(
                store, target, metric=metric, cutoff=cutoff, hops=hops,
            )
            safe_name = target.replace(" ", "_").replace("'", "")[:40]
            fname = f"neighborhood_{safe_name}.html"
            path = os.path.join(output, fname)
            fig.write_html(path)
            generated_files.append((fname, f'Neighborhood: {target}'))
            LOGGER.INFO(f"  -> {path}")

            # Threshold sweep
            LOGGER.INFO(f'Generating threshold sweep for "{target}"...')
            thresholds = None
            if thresholds_str:
                thresholds = [float(x.strip()) for x in thresholds_str.split(",")]
            fig = plot_threshold_sweep(
                store, target, metric=metric, cutoff=cutoff,
                thresholds=thresholds, hops=hops,
            )
            fname = f"threshold_sweep_{safe_name}.html"
            path = os.path.join(output, fname)
            fig.write_html(path)
            generated_files.append((fname, f'Threshold Sweep: {target}'))
            LOGGER.INFO(f"  -> {path}")

        # Index page
        _generate_index_html(generated_files, output)
        LOGGER.INFO(f"  -> {os.path.join(output, 'index.html')}")
        LOGGER.INFO(f"Generated {len(generated_files)} visualizations in {output}/")

    except KeyError as e:
        LOGGER.ERROR(str(e))
    except ImportError as e:
        LOGGER.ERROR(str(e))
    finally:
        store.close()
