import math
import itertools
import streamlit as st
import sys
from collections import OrderedDict

# Increase recursion limit for the deep optimization search
sys.setrecursionlimit(2000)

# --- CORE LOGIC FUNCTIONS (REUSED) ---

def generate_lane_splits(n_versions: int, max_total_lanes: int):
    """Generates all valid ways to assign lanes (>= 1) to n_versions."""
    if n_versions == 1:
        for l in range(1, max_total_lanes + 1):
            yield (l,)
        return

    max_l1 = max_total_lanes - (n_versions - 1)
    
    for l1 in range(1, max_l1 + 1):
        for remaining_split in generate_lane_splits(n_versions - 1, max_total_lanes - l1):
            yield (l1,) + remaining_split


def calculate_gang_waste(gang_plan: list[tuple[str, int, int]], die_repeat: float, max_lanes: int) -> float:
    """Calculates the Total Wasted Footage for a single gang, including Running and Setup Waste."""
    if not gang_plan:
        return 0.0

    # 1. Total Run Repeats (R)
    run_repeats = 0
    for _, qty, lanes in gang_plan:
        if lanes > 0:
            run_repeats = max(run_repeats, math.ceil(qty / lanes))

    # 2. Total Wasted Labels (WL)
    total_wasted_labels = 0
    total_lanes_used = sum(lanes for _, _, lanes in gang_plan)

    for _, qty, lanes in gang_plan:
        if lanes > 0:
            total_wasted_labels += (run_repeats * lanes) - qty

    # 3. Wasted Repeats (WR)
    if total_lanes_used == 0: return float('inf') 
    wasted_repeats = total_wasted_labels / total_lanes_used

    # 4a. Running Waste (WF_running)
    wf_running = (wasted_repeats * die_repeat) / 12

    # 4b. Setup Waste (WF_setup) - 50 * Die Repeat / 12 added per print job
    wf_setup = (50 * die_repeat) / 12
    
    # 4c. Total Gang Waste
    wasted_footage = wf_running + wf_setup

    return wasted_footage


def find_optimal_ganging(
    required_quantities: dict[str, int], 
    max_lanes: int, 
    die_repeat: float
) -> tuple[float, list]:
    
    memo = {}

    def solve(current_versions: tuple) -> tuple[float, list]:
        if not current_versions:
            return 0.0, []
        if current_versions in memo:
            return memo[current_versions]

        min_waste = float('inf')
        best_overall_plan = []

        for i in range(1, len(current_versions) + 1):
            for candidate_gang_subset in itertools.combinations(current_versions, i):
                
                num_versions_in_gang = len(candidate_gang_subset)
                
                if num_versions_in_gang > max_lanes:
                    continue
                
                best_gang_waste_for_subset = float('inf')
                best_lane_split = None

                for lane_split in generate_lane_splits(num_versions_in_gang, max_lanes):
                    
                    gang_plan_candidate = [(v_id, v_qty, lane_split[idx]) 
                                           for idx, (v_id, v_qty) in enumerate(candidate_gang_subset)]
                    
                    gang_waste = calculate_gang_waste(gang_plan_candidate, die_repeat, max_lanes)

                    if gang_waste < best_gang_waste_for_subset:
                        best_gang_waste_for_subset = gang_waste
                        best_lane_split = lane_split

                if best_lane_split is None:
                    continue
                
                gang_version_ids = {v_id for v_id, _ in candidate_gang_subset}
                
                remaining_requirements = tuple(
                    (v_id, v_qty) 
                    for v_id, v_qty in current_versions 
                    if v_id not in gang_version_ids
                )
                
                remaining_waste, remaining_plan = solve(remaining_requirements)

                total_current_waste = best_gang_waste_for_subset + remaining_waste

                if total_current_waste < min_waste:
                    min_waste = total_current_waste
                    
                    current_gang_info = {
                        'versions': [v_id for v_id, _ in candidate_gang_subset],
                        'lanes': list(best_lane_split),
                        'waste': best_gang_waste_for_subset
                    }
                    best_overall_plan = [current_gang_info] + remaining_plan
        
        memo[current_versions] = (min_waste, best_overall_plan)
        return min_waste, best_overall_plan


    initial_requirements_tuple = tuple(sorted(required_quantities.items()))

    final_min_waste, final_best_plan = solve(initial_requirements_tuple)
    
    return final_min_waste, final_best_plan

# --- STREAMLIT MAIN FUNCTION ---

def main():
    """The main Streamlit application interface."""
    st.set_page_config(page_title="Ganging Optimizer", layout="wide")
    st.title("Ganging Optimization Tool ⚙️")

    # --- INPUT SIDEBAR ---
    st.sidebar.header("1. Tooling Parameters")
    
    # Tooling Constants
    max_lanes = st.sidebar.number_input(
        "Max Lanes Across (Total available lanes)", 
        min_value=1, value=8, step=1, key="max_lanes"
    )
    die_repeat = st.sidebar.number_input(
        "Die Repeat (in inches)", 
        min_value=0.01, value=4.95, step=0.01, key="die_repeat"
    )

    # --- ORDER INPUT ---
    st.header("2. Enter Customer Orders")
    input_text = st.text_area(
        "Enter Version Quantities (One entry per line: Version_Name, Quantity)",
        value="V1, 10000\nV2, 300000\nV3, 5000\nV4, 15000",
        height=200,
        help="Example format: V_A, 10000"
    )
    
    # --- EXECUTION ---
    if st.button("Calculate Optimal Plan", type="primary"):
        required_quantities = OrderedDict()
        
        # 3. Parsing Logic
        try:
            lines = [line.strip() for line in input_text.split('\n') if line.strip()]
            for line in lines:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) != 2:
                    st.error("Input format error: Each line must contain exactly one comma separating the name and quantity.")
                    return

                version_name = parts[0]
                quantity = int(parts[1])
                
                if quantity <= 0:
                    st.error(f"Quantity for {version_name} must be positive.")
                    return
                
                required_quantities[version_name] = quantity
                
        except ValueError:
            st.error("Input error: Quantity must be a whole number.")
            return

        # 4. Run Optimization
        if not required_quantities:
            st.warning("Please enter at least one version quantity.")
            return
            
        with st.spinner('Running optimization, please wait... (May take time for many versions)'):
            try:
                total_waste, optimal_plan = find_optimal_ganging(
                    required_quantities=required_quantities,
                    max_lanes=max_lanes,
                    die_repeat=die_repeat
                )
            except RecursionError:
                st.error("Optimization failed: Too many versions for this method. Try fewer than 7 versions.")
                return

        # 5. Display Results
        st.subheader("3. Optimal Ganging Plan")
        st.success(f"Total Wasted Footage (Minimized): **{total_waste:.2f} ft**")

        for i, gang in enumerate(optimal_plan):
            gang_plan_for_calc = [(v_id, required_quantities[v_id], l) 
for i, gang in enumerate(optimal_plan):
            versions_list = ', '.join([f"{v} ({l} lanes)" for v, l in zip(gang['versions'], gang['lanes'])])
            
            # Recalculate Run Repeats for display
            gang_plan_for_calc = [(v_id, required_quantities[v_id], l) 
                                  for v_id, l in zip(gang['versions'], gang['lanes'])]
            
            run_repeats = max(math.ceil(qty / lanes) for _, qty, lanes in gang_plan_for_calc if lanes > 0)
            
            setup_waste = (50 * die_repeat) / 12
            running_waste = gang['waste'] - setup_waste

            st.markdown(f"**Gang {i+1}**")
            st.info(f"Run Repeats (Impressions): **{run_repeats}**")
            
            # Display versions and lanes in a table
            data = {'Version': gang['versions'], 'Lanes': gang['lanes']}
            st.table(data)
            
            # --- START FIX HERE ---
            st.markdown(f"Running Waste (from labels): **{running_waste:.2f} ft**")
            st.markdown(f"Setup Waste (50 * DR): **{setup_waste:.2f} ft**")
            st.markdown(f"**Total Gang Cost: {gang['waste']:.2f} ft**")
            # --- END FIX HERE ---
            
            st.markdown("---")


if __name__ == "__main__":

    main()
