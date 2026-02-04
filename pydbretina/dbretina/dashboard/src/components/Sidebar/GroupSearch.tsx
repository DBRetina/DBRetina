import { useState, useRef } from "react";
import { useDashboard } from "../../state/context";
import { searchGroups } from "../../api";
import { GroupMatch } from "../../state/types";

export default function GroupSearch() {
  const { state, dispatch } = useDashboard();
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<GroupMatch[]>([]);
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
        const res = await searchGroups(val, 10);
        setResults(res.matches);
        setShowDropdown(res.matches.length > 0);
      } catch {
        setResults([]);
      }
    }, 250);
  }

  function handleSelect(name: string) {
    setQuery(name);
    setShowDropdown(false);
    dispatch({ type: "SET_FOCUS_GROUP", group: name });
  }

  function handleClear() {
    setQuery("");
    setResults([]);
    setShowDropdown(false);
    dispatch({ type: "SET_FOCUS_GROUP", group: null });
  }

  return (
    <div className="sidebar-section">
      <label>
        Focus Group
        {state.focusGroup && (
          <span
            onClick={handleClear}
            style={{ cursor: "pointer", color: "var(--danger)", marginLeft: 8, fontSize: 10, textTransform: "none" }}
          >
            clear
          </span>
        )}
      </label>
      <div className="search-container">
        <input
          type="text"
          placeholder="Search diseases..."
          value={query}
          onChange={(e) => handleInput(e.target.value)}
          onFocus={() => results.length > 0 && setShowDropdown(true)}
          onBlur={() => setTimeout(() => setShowDropdown(false), 200)}
          aria-label="Search for a group by name"
          aria-autocomplete="list"
          aria-controls="group-search-results"
          role="combobox"
          aria-expanded={results.length > 0}
        />
        {showDropdown && (
          <div
            className="search-dropdown"
            id="group-search-results"
            role="listbox"
            aria-label="Search results"
          >
            {results.map((m) => (
              <div
                key={m.id}
                className="search-dropdown-item"
                onMouseDown={() => handleSelect(m.name)}
              >
                {m.name}
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
