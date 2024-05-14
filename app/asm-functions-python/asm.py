def find_asm(
    x,
    wilcoxon_to_use,
    max_p_value,
    min_nb_cpg_same_direction,
    min_nb_consecutive_cpg_same_direction,
    min_asm_region_effect,
):
    if (
        x[wilcoxon_to_use] <= max_p_value
        and abs(x["read_asm_effect"]) >= min_asm_region_effect
        and x["consecutive_sig_cpgs"] >= min_nb_consecutive_cpg_same_direction
        and x["total_sig_cpgs"] >= min_nb_cpg_same_direction
    ):
        return 1
    return 0


def find_max_consecutive_positions(qualifying_cpg_pos, all_cpg_pos):
    qualifying_lookup = {value: i for i, value in enumerate(qualifying_cpg_pos)}
    # To track the maximum length of consecutive elements
    max_consecutive = 0
    current_streak = 0
    # Iterate over all_cpg_pos to count consecutive qualifying elements
    for pos in all_cpg_pos:
        if pos in qualifying_lookup:
            if current_streak == 0:  # Starting a new consecutive sequence
                current_streak = 1
            else:
                # Check if this is truly consecutive in terms of qualifying_cpg_pos
                last_pos = all_cpg_pos.index(pos) - 1
                if last_pos >= 0 and all_cpg_pos[last_pos] in qualifying_lookup:
                    current_streak += 1
                else:
                    current_streak = 1
        else:
            if current_streak > max_consecutive:
                max_consecutive = current_streak
            current_streak = 0
    # Final check to update max_consecutive at the end of the loop
    if current_streak > max_consecutive:
        max_consecutive = current_streak
    return max_consecutive
