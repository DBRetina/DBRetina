import { useState, useRef } from "react";
import { useDashboard } from "../../state/context";
import { searchFeatures } from "../../api";

export default function FeatureSearch() {
  const { dispatch } = useDashboard();
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<string[]>([]);
  const [showDropdown, setShowDropdown] = useState(false);
  const debounceRef = useRef<ReturnType<typeof setTimeout>>();

  function handleInput(val: string) {
    setQuery(val);
    clearTimeout(debounceRef.current);
    if (val.length < 2) {
      setResults([]);
      setShowDropdown(false);
      return;
    }
    debounceRef.current = setTimeout(async () => {
      try {
        const res = await searchFeatures(val, 10);
        setResults(res.groups);
        setShowDropdown(res.groups.length > 0);
      } catch {
        setResults([]);
      }
    }, 250);
  }

  function handleSelect(groupName: string) {
    setShowDropdown(false);
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
          onFocus={() => results.length > 0 && setShowDropdown(true)}
          onBlur={() => setTimeout(() => setShowDropdown(false), 200)}
        />
        {showDropdown && (
          <div className="search-dropdown">
            {results.map((name) => (
              <div
                key={name}
                className="search-dropdown-item"
                onMouseDown={() => handleSelect(name)}
              >
                {name}
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
