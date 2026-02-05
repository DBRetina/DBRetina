import { useState, useRef } from "react";
import { useDashboard } from "../../state/context";
import { searchGenes, searchFeatures } from "../../api";

interface GeneMatch {
  gene: string;
  group_count: number;
}

export default function FeatureSearch() {
  const { dispatch } = useDashboard();
  const [query, setQuery] = useState("");
  const [geneResults, setGeneResults] = useState<GeneMatch[]>([]);
  const [groupResults, setGroupResults] = useState<string[]>([]);
  const [selectedGene, setSelectedGene] = useState<string | null>(null);
  const [showDropdown, setShowDropdown] = useState(false);
  const debounceRef = useRef<ReturnType<typeof setTimeout>>();

  function handleInput(val: string) {
    setQuery(val);
    setSelectedGene(null);
    setGroupResults([]);
    clearTimeout(debounceRef.current);
    if (val.length < 2) {
      setGeneResults([]);
      setShowDropdown(false);
      return;
    }
    debounceRef.current = setTimeout(async () => {
      try {
        const res = await searchGenes(val, 15);
        setGeneResults(res.genes);
        setShowDropdown(res.genes.length > 0);
      } catch {
        setGeneResults([]);
      }
    }, 200);
  }

  async function handleSelectGene(gene: string) {
    setQuery(gene);
    setSelectedGene(gene);
    setGeneResults([]);
    setShowDropdown(false);

    // Now fetch groups containing this gene
    try {
      const res = await searchFeatures(gene, 20);
      setGroupResults(res.groups);
    } catch {
      setGroupResults([]);
    }
  }

  function handleSelectGroup(groupName: string) {
    dispatch({ type: "SET_FOCUS_GROUP", group: groupName });
  }

  return (
    <div className="sidebar-section">
      <label>Gene Search</label>
      <div className="search-container">
        <input
          type="text"
          placeholder="Search by gene (e.g. TREM2)..."
          value={query}
          onChange={(e) => handleInput(e.target.value)}
          onFocus={() => geneResults.length > 0 && setShowDropdown(true)}
          onBlur={() => setTimeout(() => setShowDropdown(false), 200)}
        />
        {showDropdown && geneResults.length > 0 && (
          <div className="search-dropdown">
            <div
              style={{
                padding: "4px 8px",
                fontSize: 10,
                color: "var(--text-secondary)",
                borderBottom: "1px solid var(--border-color)",
              }}
            >
              Genes matching &ldquo;{query}&rdquo;
            </div>
            {geneResults.map((match) => (
              <div
                key={match.gene}
                className="search-dropdown-item"
                onMouseDown={() => handleSelectGene(match.gene)}
                style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}
              >
                <span>{match.gene}</span>
                <span style={{ fontSize: 10, color: "var(--text-secondary)" }}>
                  {match.group_count} groups
                </span>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Groups containing the selected gene */}
      {selectedGene && groupResults.length > 0 && (
        <div style={{ marginTop: 8 }}>
          <div
            style={{
              fontSize: 10,
              color: "var(--text-secondary)",
              marginBottom: 4,
            }}
          >
            Groups containing <strong>{selectedGene}</strong> ({groupResults.length})
          </div>
          <div
            style={{
              maxHeight: 160,
              overflowY: "auto",
              border: "1px solid var(--border-color)",
              borderRadius: 4,
            }}
          >
            {groupResults.map((name, i) => (
              <div
                key={name}
                className="search-dropdown-item"
                onClick={() => handleSelectGroup(name)}
                style={{
                  borderBottom:
                    i < groupResults.length - 1
                      ? "1px solid var(--border-color)"
                      : "none",
                }}
              >
                {name}
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
