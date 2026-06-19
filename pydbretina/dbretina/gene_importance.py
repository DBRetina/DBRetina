"""Gene-level importance scoring for DBRetina pairwise data."""

import json
import math
from collections import defaultdict
from typing import Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .pairwise_store import PairwiseStore


class GeneImportance:
    """Score individual genes by their importance in driving disease similarity.

    Supports three scoring methods:
    - **edge_weighted**: Sum of edge weights (metric values) per gene.
    - **hypergraph**: Local frequency x inverse global frequency (TF-IDF style).
    - **projection**: PageRank on gene co-occurrence graph.

    Usage::

        store = PairwiseStore("/path/to/pairwise", dbri_path="/path/to/index.dbri")
        gi = GeneImportance(store)

        # Explain why two diseases are similar
        df = gi.explain_pair("autism spectrum disorders", "schizophrenia")

        # Find hub genes for a disease
        hubs = gi.hub_genes("autism spectrum disorders", hops=2)

        # Save results
        scores = gi.projection_scores("autism spectrum disorders", hops=2)
        gi.to_parquet("gene_scores.parquet", scores,
                      {"target": "autism spectrum disorders", "method": "projection"})
    """

    def __init__(self, store: PairwiseStore, graph=None):
        """Create a GeneImportance scorer.

        Args:
            store: PairwiseStore with dbri_path set (required for gene access).
            graph: Optional pre-built PairwiseGraph for neighborhood extraction.
                   If None, built lazily when hops-based methods are called.
        """
        self._store = store
        self._graph = graph

    def _get_graph(self):
        """Get or build the PairwiseGraph for neighborhood extraction."""
        if self._graph is not None:
            return self._graph
        from .pairwise_graph import PairwiseGraph
        self._graph = PairwiseGraph(self._store, metric="ochiai", cutoff=10)
        return self._graph

    def _get_neighborhood_diseases(self, group_name: str, hops: int) -> list[str]:
        """Get disease names in the N-hop neighborhood of a group."""
        graph = self._get_graph()
        sub = graph.subgraph_around(group_name, hops=hops)
        names = list(sub._names_map.values())
        sub.close()
        return names

    # ── Scoring Methods ────────────────────────────────────────────

    def edge_weighted_scores(
        self, group_name: str, metric: str = "ochiai", cutoff: float = 20,
        normalize: bool = True
    ) -> dict[str, float]:
        """Score genes by sum of edge weights they participate in.

        For each pair involving group_name where metric >= cutoff,
        each shared gene accumulates the metric value of that pair.
        When normalize=True (default), scores are divided by gene degree
        to prevent hub genes from dominating.

        Args:
            group_name: Target disease/group name.
            metric: Similarity metric to use as weight.
            cutoff: Minimum metric value to include a pair.
            normalize: If True, normalize by gene degree (default: True).

        Returns:
            Dict mapping gene_name -> accumulated score.
        """
        if metric not in PairwiseStore.METRICS:
            raise ValueError(
                f"Unknown metric '{metric}'. Valid: {', '.join(PairwiseStore.METRICS)}"
            )
        store = self._store
        gene_sets = store._load_gene_sets()
        target_key = group_name.lower()

        if target_key not in gene_sets:
            available = list(gene_sets.keys())[:5]
            raise KeyError(
                f"Group '{group_name}' not found in gene sets. "
                f"Available groups include: {', '.join(available)}{'...' if len(gene_sets) > 5 else ''}"
            )

        target_genes = set(gene_sets[target_key])
        gid = store.group_id(group_name)
        if gid is None:
            raise KeyError(f"Group not found: {group_name}")

        # Get all pairs involving this group above the cutoff
        df = store.sql(
            f"SELECT group_1_id, group_2_id, {metric} "
            f"FROM pairs "
            f"WHERE (group_1_id = {gid} OR group_2_id = {gid}) "
            f"AND {metric} >= {cutoff}"
        ).fetchdf()

        names_map = store.get_names_map()
        raw_scores: dict[str, float] = defaultdict(float)
        gene_degree: dict[str, int] = defaultdict(int)

        for _, row in df.iterrows():
            other_id = row["group_2_id"] if row["group_1_id"] == gid else row["group_1_id"]
            other_name = names_map.get(other_id, "").lower()
            if other_name not in gene_sets:
                continue

            other_genes = set(gene_sets[other_name])
            shared = target_genes & other_genes
            weight = float(row[metric])

            for gene in shared:
                raw_scores[gene] += weight
                gene_degree[gene] += 1

        # Normalize by gene degree to prevent hub gene dominance
        if normalize:
            scores = {
                gene: raw_scores[gene] / gene_degree[gene]
                for gene in raw_scores
                if gene_degree[gene] > 0
            }
        else:
            scores = dict(raw_scores)

        return scores

    def hypergraph_scores(
        self, group_name: str, hops: int = 2
    ) -> dict[str, float]:
        """Score genes by local frequency x inverse global frequency (TF-IDF).

        For each gene in the target disease:
          TF = local_count / neighborhood_size (term frequency)
          IDF = log(total_groups / global_count) (inverse document frequency)
          score = TF × IDF (standard TF-IDF formula)

        Genes specific to the disease neighborhood score high.
        Ubiquitous genes (TNF, IL6) get down-weighted by IDF.

        Args:
            group_name: Target disease/group name.
            hops: Number of hops for neighborhood extraction.

        Returns:
            Dict mapping gene_name -> TF-IDF score.
        """
        store = self._store
        gene_sets = store._load_gene_sets()
        gene_index = store._get_gene_index()
        total_groups = store.num_groups
        target_key = group_name.lower()

        if target_key not in gene_sets:
            available = list(gene_sets.keys())[:5]
            raise KeyError(
                f"Group '{group_name}' not found in gene sets. "
                f"Available groups include: {', '.join(available)}{'...' if len(gene_sets) > 5 else ''}"
            )

        target_genes = set(gene_sets[target_key])

        # Get neighborhood diseases
        neighborhood = self._get_neighborhood_diseases(group_name, hops)

        # Count how many neighborhood diseases contain each target gene
        neighborhood_sets = {}
        for disease in neighborhood:
            dk = disease.lower()  # Fixed: simplified case handling
            if dk in gene_sets:
                neighborhood_sets[dk] = set(gene_sets[dk])

        n_neighborhood = len(neighborhood_sets)
        scores: dict[str, float] = {}
        for gene in target_genes:
            local_count = sum(
                1 for gs in neighborhood_sets.values() if gene in gs
            )
            global_count = len(gene_index.get(gene, set()))
            if global_count > 0 and local_count > 0:
                # TF: fraction of neighborhood diseases containing gene
                tf = local_count / n_neighborhood if n_neighborhood > 0 else 0
                # IDF: penalize genes appearing in many diseases globally
                # Standard TF-IDF formula: log(N/df) where N=total docs, df=doc frequency
                idf = math.log(total_groups / global_count)
                # Final score: standard TF × IDF
                # (Previous version had non-standard enrichment × IDF which was
                # equivalent to TF × (N/df) × log(N/df), creating quadratic penalty)
                scores[gene] = tf * idf

        return scores

    def projection_scores(
        self, group_name: str, hops: int = 2
    ) -> dict[str, float]:
        """Score genes via PageRank on gene co-occurrence graph.

        Builds a gene-gene graph from the N-hop neighborhood where
        two genes are connected if they co-occur in the same disease.
        Edge weight = number of diseases sharing both genes.
        Runs PageRank on this graph.

        Args:
            group_name: Target disease/group name.
            hops: Number of hops for neighborhood extraction.

        Returns:
            Dict mapping gene_name -> PageRank score.
        """
        try:
            import igraph as ig
        except ImportError:
            raise ImportError(
                "igraph is required for projection_scores(). "
                "Install with: pip install igraph or pip install 'DBRetina[all]'"
            )

        store = self._store
        gene_sets = store._load_gene_sets()
        target_key = group_name.lower()

        if target_key not in gene_sets:
            available = list(gene_sets.keys())[:5]
            raise KeyError(
                f"Group '{group_name}' not found in gene sets. "
                f"Available groups include: {', '.join(available)}{'...' if len(gene_sets) > 5 else ''}"
            )

        target_genes = set(gene_sets[target_key])

        # Get neighborhood diseases
        neighborhood = self._get_neighborhood_diseases(group_name, hops)

        # Collect gene sets for neighborhood diseases
        neighborhood_gene_sets = []
        for disease in neighborhood:
            dk = disease.lower()  # Fixed: simplified case handling
            if dk in gene_sets:
                # Only include genes that are also in the target disease
                shared = target_genes & set(gene_sets[dk])
                if shared:
                    neighborhood_gene_sets.append(shared)

        if not neighborhood_gene_sets:
            return {}

        # Build gene list and index
        all_genes = sorted(target_genes)
        gene_idx = {g: i for i, g in enumerate(all_genes)}

        # Build gene-gene co-occurrence graph
        edge_weights: dict[tuple[int, int], int] = defaultdict(int)
        for gene_set in neighborhood_gene_sets:
            genes_in_set = sorted(gene_set)
            for i in range(len(genes_in_set)):
                gi = gene_idx.get(genes_in_set[i])
                if gi is None:
                    continue
                for j in range(i + 1, len(genes_in_set)):
                    gj = gene_idx.get(genes_in_set[j])
                    if gj is None:
                        continue
                    edge_weights[(gi, gj)] += 1

        # Build igraph
        g = ig.Graph(n=len(all_genes), directed=False)
        edges = list(edge_weights.keys())
        weights = [edge_weights[e] for e in edges]

        if not edges:
            # No co-occurrence edges; return uniform scores
            uniform = 1.0 / len(all_genes) if all_genes else 0
            return {gene: uniform for gene in all_genes}

        g.add_edges(edges)
        g.es["weight"] = weights

        pr = g.pagerank(damping=0.85, weights="weight")
        return {all_genes[i]: pr[i] for i in range(len(all_genes))}

    # ── Convenience Methods ────────────────────────────────────────

    def explain_pair(
        self, group_a: str, group_b: str, method: str = "hypergraph",
        metric: str = "ochiai", cutoff: float = 20,
    ) -> pd.DataFrame:
        """Rank shared genes between two diseases by importance.

        Args:
            group_a: First disease/group name.
            group_b: Second disease/group name.
            method: Scoring method ("hypergraph", "edge_weighted", or "projection").
            metric: Similarity metric for the edge_weighted method (ignored otherwise).
            cutoff: Minimum metric value for the edge_weighted method (ignored otherwise).

        Returns:
            DataFrame with columns: gene, score, local_freq, global_freq, specificity.
        """
        store = self._store
        shared = store.shared_features(group_a, group_b)
        gene_index = store._get_gene_index()
        total_groups = store.num_groups

        # Get scores from selected method for context
        if method == "edge_weighted":
            all_scores = self.edge_weighted_scores(group_a, metric=metric, cutoff=cutoff)
        elif method == "projection":
            all_scores = self.projection_scores(group_a)
        else:
            all_scores = self.hypergraph_scores(group_a)

        rows = []
        for gene in shared:
            global_freq = len(gene_index.get(gene, set()))
            specificity = 1.0 / global_freq if global_freq > 0 else 0
            score = all_scores.get(gene, 0.0)
            rows.append({
                "gene": gene,
                "score": score,
                "global_freq": global_freq,
                "specificity": round(specificity, 6),
            })

        df = pd.DataFrame(rows)
        if len(df) > 0:
            df = df.sort_values("score", ascending=False).reset_index(drop=True)
        return df

    def hub_genes(
        self,
        group_name: str,
        hops: int = 2,
        method: str = "projection",
        top_n: int = 30,
        metric: str = "ochiai",
        cutoff: float = 20,
    ) -> pd.DataFrame:
        """Rank all genes in a disease's neighborhood by importance.

        Args:
            group_name: Target disease/group name.
            hops: Number of hops for neighborhood.
            method: Scoring method ("projection", "hypergraph", or "edge_weighted").
            top_n: Number of top genes to return.
            metric: Similarity metric for the edge_weighted method (ignored otherwise).
            cutoff: Minimum metric value for the edge_weighted method (ignored otherwise).

        Returns:
            DataFrame with columns: gene, score, num_diseases.
        """
        if method == "edge_weighted":
            scores = self.edge_weighted_scores(group_name, metric=metric, cutoff=cutoff)
        elif method == "hypergraph":
            scores = self.hypergraph_scores(group_name, hops=hops)
        else:
            scores = self.projection_scores(group_name, hops=hops)

        gene_index = self._store._get_gene_index()
        rows = []
        for gene, score in sorted(scores.items(), key=lambda x: -x[1])[:top_n]:
            rows.append({
                "gene": gene,
                "score": round(score, 6),
                "num_diseases": len(gene_index.get(gene, set())),
            })
        return pd.DataFrame(rows)

    # ── Persistence ────────────────────────────────────────────────

    @staticmethod
    def to_parquet(
        path: str, scores: dict[str, float], metadata: Optional[dict] = None
    ) -> None:
        """Save gene importance scores to a Parquet file.

        Args:
            path: Output file path.
            scores: Dict mapping gene_name -> score.
            metadata: Optional metadata dict stored in Parquet schema metadata.
        """
        df = pd.DataFrame([
            {"gene": g, "score": s} for g, s in scores.items()
        ]).sort_values("score", ascending=False).reset_index(drop=True)

        table = pa.Table.from_pandas(df)
        if metadata:
            existing = table.schema.metadata or {}
            existing[b"dbretina_gene_importance"] = json.dumps(metadata).encode()
            table = table.replace_schema_metadata(existing)
        pq.write_table(table, path)

    @staticmethod
    def from_parquet(path: str) -> pd.DataFrame:
        """Load gene importance scores from a Parquet file.

        Args:
            path: Path to a previously saved scores file.

        Returns:
            DataFrame with gene and score columns.
        """
        return pq.read_table(path).to_pandas()

    def __repr__(self):
        graph_info = f", graph={self._graph}" if self._graph else ""
        return f"GeneImportance(store={self._store}{graph_info})"
