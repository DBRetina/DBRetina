import { useCallback } from "react";
import { useDashboard } from "../../state/context";
import DatasetInfo from "./DatasetInfo";
import MetricSelector from "./MetricSelector";
import CutoffSlider from "./CutoffSlider";
import GroupSearch from "./GroupSearch";
import FeatureSearch from "./FeatureSearch";
import GeneExplorer from "./GeneExplorer";
import HubGenesPanel from "./HubGenesPanel";
import PathFinder from "./PathFinder";
import CommunityPanel from "./CommunityPanel";
import ClusteringPanel from "./ClusteringPanel";
import ExportPanel from "./ExportPanel";
import { FilterBuilder, FilterCondition } from "./AdvancedFilter";

export default function Sidebar() {
  const { state, dispatch } = useDashboard();

  // Handle filter application (placeholder - can be extended to update graph)
  const handleFilterApply = useCallback(
    (filters: FilterCondition[], logic: "AND" | "OR") => {
      console.log("Applying filters:", filters, "with logic:", logic);
      // TODO: Integrate with graph state to show filtered results
      // This could dispatch an action to update a "filters" state
    },
    []
  );

  const handleFilterClear = useCallback(() => {
    console.log("Clearing filters");
    // TODO: Clear filter state
  }, []);

  return (
    <div className="sidebar">
      <div className="sidebar-header">
        <h1>DBRetina</h1>
      </div>
      <DatasetInfo />
      <MetricSelector />
      <CutoffSlider />
      <GroupSearch />
      {state.info?.has_genes && <FeatureSearch />}
      {state.info?.has_genes && (
        <>
          <GeneExplorer />
          <HubGenesPanel />
          <PathFinder />
        </>
      )}
      <FilterBuilder onApply={handleFilterApply} onClear={handleFilterClear} />
      <ClusteringPanel />
      <ExportPanel />
      <CommunityPanel />
    </div>
  );
}
